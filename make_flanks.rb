require 'optparse'

# Convert distance strings like "50u" or "20d" to integer
# Negative for upstream, positive for downstream
def parse_distance(s)
  return 0 if s == '0'

  match = s.downcase.match(/^(\d+)([ud])$/)
  raise "distance doesn't match format: #{s}" unless match

  match[2] == 'u' ? -Integer(match[1]) : Integer(match[1])
end

# Parse command-line arguments
OptionParser.new do |opts|
  opts.banner = "Usage: make_flanks.rb <upstream> <downstream> <input.bed>"
end.parse!

raise "Usage: make_flanks.rb <upstream> <downstream> <input.bed>" unless ARGV.length == 3

flank_5 = parse_distance(ARGV[0])
flank_3 = parse_distance(ARGV[1])
input_fn = ARGV[2]

# Read input lines from file or STDIN
lines = ['stdin', '-'].include?(input_fn.downcase) ? STDIN.readlines : File.readlines(input_fn)

lines.each do |line|
  next if line.strip.empty? || line.start_with?('#')
  row = line.strip.split("\t")

  chr    = row[0]
  start  = Integer(row[1])
  stop   = Integer(row[2])
  name   = row[3]
  strand = row[5]

  # Use start as center (1-nt length)
  center = start

  from, to = case strand
  when '+'
    [center + flank_5, center + flank_3 + 1]
  when '-'
    [center - flank_3, center - flank_5 + 1]
  else
    raise "Unknown strand: #{strand}"
  end

  from = [from, 0].max  # ensure coordinates are not negative
  puts [chr, from, to, name, '.', strand].join("\t")
end
