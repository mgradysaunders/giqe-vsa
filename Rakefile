require 'rake'

desc "Default task."
task :default => :bin

desc "Build binary."
task :bin do
    mkdir_p "bin"
    sh "g++-8 -std=c++17 -Wall -Wextra -O3 -Ilib/pcg-cpp/include -DNDEBUG src/main.cpp -o bin/giqe-vsa"
end
