#!/usr/bin/env ruby 

require 'fileutils'

DIR = ARGV[0]

class Table

  attr_reader :c
  attr_reader :cnames

  def initialize( fname )
    read(fname)
  end

  def read( fname )
    @c = nil
    @cnames = {}

    File.open(fname, 'r') do |file|
       file.each do |line|
         row = line.split

		 if  @cnames.empty? then
			row.each_with_index{ |n,i| @cnames[n.to_sym] = i }
		 	next
		 end
        
         @c = Array.new(row.count){[]}  unless @c
    
         row.each_with_index { |val,k| @c[k] << val.to_f } 
       end
    end  
  end 

  def col( n )
    @c[@cnames[n]]
  end

  def cols?
    @c.count
  end

end

t = Table.new( File.join(DIR, 'stats.txt') )

max_nodes = t.col( :actNds ).max.to_i
max_parts = t.col( :nPart ).max.to_i

puts "Max active nodes: #{max_nodes}"
puts "Max particles: #{max_parts}"

n_frames = t.col( :frame )[-1]+1
n_steps  = t.col( :frame ).count
tot_time = t.col( :totTime ).reduce(:+)
slv_time = t.col( :slvTime ).reduce(:+)
asm_time = t.col( :asmTime ).reduce(:+)

puts "Tot time: #{tot_time}"
puts "Tot frames: #{n_frames}"
puts "Avg time per frame: #{tot_time/n_frames}"
puts "Avg steps per frame: #{n_steps/n_frames}"
puts "Percent solve: #{slv_time*100/tot_time}"
puts "Percent assembly: #{asm_time*100/tot_time}"


config_file = File.join(DIR,'config')
fps = 0
File.open( config_file ).each_line{ |l|
	if l =~ /^res\s+([-0-9\s.eE]+)\s*$/
          puts "Grid: #{$1}"
	end
        if l =~ /^fps\s+([-0-9.eE]+)\s*$/
		fps = $1.to_f
	end
}
puts "Target fps: #{fps}"
puts "Real time ratio: #{fps*tot_time/n_frames}"

