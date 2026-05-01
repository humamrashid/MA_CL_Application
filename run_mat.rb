#!/usr/bin/ruby
# Run script for matrix multiplier.

def run_mat(b, m, s, t, l, u)
  count = 0
  (1..m).each do |n|
    n = b**n if b != 0
    system("echo #{n} > in_mat.txt")
    system("./rand_grid #{n * n} #{l} #{u} #{n} >> in_mat.txt")
    system("./rand_grid #{n * n} #{l} #{u} #{n} >> in_mat.txt")
    puts if count > 0
    puts "*** Run ##{count += 1} (#{n} x #{n} matrix) ***\n"
    case t
    when 'a' then system("./mat_mul #{s} in_mat.txt")
    when 'e' then system("./mat_mul #{s} in_mat.txt | grep -i elapsed")
    when 'o' then system("./mat_mul #{s} in_mat.txt | grep -i number")
    end
  end
end

if ARGV.length != 6
  abort "Usage: #{$PROGRAM_NAME} -<n|s> -<a|e|o> -<#> <bound> <l> <u>"
end
s = ARGV[0]
case s
when '-n' then method = 'Naive'
when '-s' then method = 'Strassen'
else
  abort "Invalid method option: #{s}"
end
t = ARGV[1][1]
abort "Invalid output option: -#{t}" if t != "a" and t != "e" and t != "o"
b = ARGV[2][1].to_i
abort 'Power base must be >= 0' if b < 0
m = ARGV[3].to_i
abort 'Upper bound must be > 0' if m <= 0
l = ARGV[4].to_i
abort 'Lower bound for random values must be >= 0' if l < 0
u = ARGV[5].to_i
abort 'Upper bound for random values must be >= lower bound' if u < l
puts "Method: #{method}"
run_mat(b, m, s, t, l, u)

# EOF.
