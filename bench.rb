#!/usr/bin/env ruby

require 'open3'
require 'fileutils'

BENCH_DIR = ARGV[0] || 'bench'
BIN_DIR   = ARGV[1] || 'apps'

OPTIONS = [
	'-v', 'Error',
	'-o',
	'-f', 1e4,
	'-e', 0.01
]

SOLVERS = [ 
	['-a', 0, '-n', 1e6],
	['-a', 1],
	['-a', 3, '-g', 2 ],
	['-a', 3, '-g', 3 ],
	['-a', 3, '-g', 4 ],
]

TOLERANCES  = [1e-6,1e-7]
TOL_TIGHTEN = 1e-1

INCREMENTAL = true

MAX_TIME = 1e99

def each_subdir dir
  Dir.foreach( dir ) do |file|
    next if file.start_with? '.' 
    next unless Dir.exists? File.join( dir, file )
    yield file
  end
end

def solve( src_dir, out_dir, file, tol, args )
	dest = File.join( out_dir, "#{file}-#{args.join ''}.log" ) ;
	return if INCREMENTAL and File.exists? dest

	cmd = [ File.join(BIN_DIR, 'd6_solvePrimal'), 
		File.join( src_dir, file),
		'-t', tol*TOL_TIGHTEN ] << args << OPTIONS 

	puts cmd.join " "
	Open3.pipeline( cmd.join(" "), :out => dest )
end

def parse_variant( v )
	case v[1].to_i
		when 0 then 'standard'
		when 1 then 'descent'
		when 2 then 'conjugated'
		when 3 then 'apgd'
		when 4 then 'spg'
	end
end

def parse_solver( s )
	args = s.split '-'
	case args[0].to_i
		when 0 then 'gs'
		when 1 then 'cadoux-gs'
		when 2 then "cadoux-pg-#{parse_variant args[1]}"
		when 3 then "cadoux-pg-expr-#{parse_variant args[1]}"
	end
end

def read_solve_time( tol, file )
	File.open( file ).each_line do |l|
		vals = l.split
		return vals[1].to_f if( vals[2].to_f < tol )
	end

	return MAX_TIME
end

def gen_diagram( perfs )
	n = 0
	diag = { :t => [0], :n => [0] }
	perfs.each do |t|
		next if t > MAX_TIME * 1e-10 

		diag[:t] << t
		diag[:n] << n
		n += 1
		diag[:t] << t
		diag[:n] << n
	end

	diag
end

def merge_diagrams( perfs )
	
	glob = { :t => [], :n => {} }
	indexes = {}

	perfs.each{ |s, v| 
		glob[:n][s] = []
		indexes[s] = 0
	}

	t = 0

	loop do
		t = perfs.map{ |s,v| (v[:t][indexes[s]] || MAX_TIME) }.min 
		break if t == MAX_TIME

		glob[:t] << t

		perfs.each do |s,v|
			idx = indexes[s]
			indexes[s] += 1 if v[:t][idx] and v[:t][idx] <= t  	
			glob[:n][s] << v[:n][indexes[s] -1]
		end

	end

	glob
end

# I - Solve problems
each_subdir BENCH_DIR do |scene|
	puts "# Entering scene #{scene}"

	dir_path = File.join( BENCH_DIR, scene ) 
	Dir.foreach( dir_path ) do |file|
		next unless file.end_with? '.d6'
		
		TOLERANCES.each do |tol|
			out_dir = File.join( dir_path, "tol-#{tol}" )
			FileUtils.mkdir out_dir unless Dir.exists? out_dir

			SOLVERS.each do |slv|
				solve( dir_path, out_dir, file, tol, slv )
			end
		end
	end

end

# II - Analyze results
all_solve_times = []
global_perfs = []

TOLERANCES.each do |tol|
	puts "# Processing tolerance #{tol}"
	solve_times = {}
	solve_perfs = {}

	each_subdir BENCH_DIR do |scene|
		times = {}
		
		# Read solve times
		dir_path = File.join( BENCH_DIR, scene, "tol-#{tol}" ) 
		Dir.foreach( dir_path ) do |file|
			if file.strip =~ /^([-\w0-9]+)\.d6--a([-.\w0-9]+).log$/ then
				s = parse_solver($2)
				t = read_solve_time( tol, File.join(dir_path, file) )
				times[$1]  ||= {} 
				times[$1][s] = t
			end
		end

		next if times.count == 0
		
		# Compile perf ratios
		perfs = {}
		times.each do |f,v| 
			best_perf = v.map{|s,t| t}.min 
			v.each do |s,t| 
				perfs[s] ||= [] 
				perfs[s] << (t-best_perf)/best_perf
			end
		end
		

		solve_times[scene] = times
		solve_perfs[scene] = perfs
	end
	
	# Gather perfs across scenes 
	tol_perfs = {}
	solve_perfs.each do |scene,v| 
		v.each{ |slv, perfs|
			( tol_perfs[slv] ||= [] ).concat perfs
		}
	end
	tol_perfs.each{ |k,v| v.sort! }
	
	global_perfs << tol_perfs
	all_solve_times << solve_times
end

#puts all_solve_times
#puts global_perfs

# III - Output perf diagrams
TOLERANCES.each_with_index do |tol, idx|
	diags = {}
	global_perfs[idx].each do |slv, perfs|
		diag = gen_diagram perfs

		File.open( File.join(BENCH_DIR, "perf-#{tol}-#{slv}" ), "w" ) { |f|
			diag[:t].each_with_index{ |t, i|
				f.puts "#{t}\t#{diag[:n][i]}"
			}
		}

		diags[slv] = diag
	end

	global_diag = merge_diagrams( diags )
	File.open( File.join(BENCH_DIR, "perf-#{tol}-all" ), "w" ) { |f|
		f.puts "#t\t#{global_diag[:n].map{ |s,v| s }.join("\t")}"
		global_diag[:t].each_with_index{ |t, i|
			f.puts "#{t}\t#{global_diag[:n].map{ |s,v| global_diag[:n][s][i] }.join("\t")}"
		}
	}
end

