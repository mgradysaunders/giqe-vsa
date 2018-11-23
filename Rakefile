require 'rake'

desc "Default task."
task :default => :bin

desc "Build binary."
task :bin do
    sh "g++-8 -std=c++17 -O3 -Ilib/pcg-cpp/include -DNDEBUG src/main.cpp -o bin/giqe-vsa"
end
