#!/usr/bin/env ruby

path = ''
dest = 'CMakeSources.txt'

if ARGV.count > 0 then
  path = ARGV[0] 
  dest = File.join(path, dest)
end

#puts "git ls-files #{path} | grep -e '\\.\\(cc\\|hh\\|c\\|h\\|hpp\\|cpp\\)\\$'  > #{dest}"
`git ls-files #{path} | grep -e "\\.\\(cc\\|hh\\|c\\|h\\|hpp\\|cpp\\)\\$"  > #{dest} `

