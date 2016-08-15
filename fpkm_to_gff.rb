input = File.open("simRef_mau12Female_genes.fpkm_tracking","r")
output = File.open("fpkmFemale.gff","a")

text = input.read
lines = text.split("\n")

for line in lines[1..-1]
  output_line = []

  cells = line.split("\t")
  
  locus = cells[6]
  prefx = locus.split(":")[0]
  #chrom = prefx.split('_')[1]
  range = locus.split(':')[1]

  if range
    range_split = range.split('-')
    range_lower = range_split[0]
    range_upper = range_split[1]
  else
    range_lower = ''
    range_upper = ''
  end

  output_line << prefx
  output_line << 'ruby'
  output_line << 'gene'
  output_line << range_lower
  output_line << range_upper
  output_line << '.'
  output_line << '+'
  output_line << '.'
  output_line << cells[0]

  output_string = output_line.join("\t")
  #puts output_string
  output << output_string+"\n"
end