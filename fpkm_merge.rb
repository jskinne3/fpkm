# This script takes two fpkm_tracking files, for male and female fruit files,
# and attempts to guess when an RNA read in one sex represents the same gene
# the other sex. More than one male read may be matched with the same female.
# The guesses are output as a .csv file.

require 'csv'

PATH = '/Users/john/Desktop/bioinformatic/merging/automerge/' # test

def open(s)
  # Open an fpkm_tracking file. Function is run twice, with
  # "M" and "F" as arguments, to open two files for comparison.
  file = File.read PATH+"mau12#{s}_mauPBref_genes.fpkm_tracking"
  rnas = CSV.parse(file, headers: true, col_sep: "\t")
  #rnas = rnas.drop(rnas.length-1000) # Subset of the lines, for testing
  data = rnas.map{|r| { # Take only needed data columns
    contig:  r['locus'].split('|')[0],
    locus: r['locus'].split(':')[1],
    fpkm: r['FPKM'],
    id: r['gene_id']
  }}
  return data
end

def overlap(a, b) # Calculate percent (decimal form) overlap of two RNAs
  # Split up the locus strings into endpoints
  a0, a1 = a.split('-').map{|n| n.to_i}
  b0, b1 = b.split('-').map{|n| n.to_i}
  # Discard if two RNAs have no overlap
  (return 0.0) if ((a1 < b0) || (b1 < a0))
  # Order the four endpoints of the two RNAs by locus
  endpts = [a0, a1, b0, b1].sort!
  # Return the ratio of the inner pair / outer pair of endpoints
  return ((endpts[1] - endpts[2]).abs.to_f / (endpts[0] - endpts[3]).abs.to_f)
end

def match(item, list)
  # Function used to find one RNA from a list of RNAs
  # with the highest overlap with the RNA item to be matched.
  rna_found, percent_max = nil, 0.0
  for line in list
    if (line[:contig] == item[:contig])
      percent_found = overlap(line[:locus], item[:locus])
      # 10% minimum overlap cutoff
      if percent_found > 0.1
        # If the overlap on this pass is higher than any overlap found previously...
        if (percent_found > percent_max)
          rna_found = line # ...then the current RNA is the best match so far... 
          percent_max = percent_found # ...and the max overlap found should be increased.
        end
      end
    end
  end
  return rna_found # return an RNA or nil
end

def main
  output, count = [], 0
  m_reads, f_reads = open('M'), open('F') # Open two fpkm_tracking files
  f_copy = f_reads.dup # Copy will keep track of unmatched female leads
  # Iterate over each male RNA read and try to find the corresponding female read
  m_reads.each do |m_read|
    f_match = match(m_read, f_reads) # Returns a matching female RNA
    if f_match
      # Remove the match from the copy list so it's not added at the end
      f_copy.delete(f_match)
      output << {m: m_read, f: f_match}
      count += 1
    else
      output << {m: m_read, f: nil}
    end
  end
  puts "Matched #{count} genes"
  # Add unmatched female RNA reads to the output also
  f_copy.each do |f_read|
    output << {m: nil, f: f_read}
  end
  # Sort output rows by contig, then by locus numbers
  output = output.sort_by do |r| 
    [
      (r[:m] ? r[:m][:contig].gsub('Segkk','').to_i : r[:f][:contig].gsub('Segkk','').to_i),
      (r[:m] ? r[:m][:locus].split('-')[0].to_i : r[:f][:locus].split('-')[0].to_i)
    ]
  end
  # Format output columns
  lines = output.map do |row|
    [
      row[:m] ? row[:m][:contig] : row[:f][:contig],
      row[:m] ? row[:m][:id] : '',
      row[:m] ? row[:m][:locus] : '',
      row[:m] ? row[:m][:fpkm] : '',
      row[:f] ? row[:f][:id] : '',
      row[:f] ? row[:f][:locus] : '',
      row[:f] ? row[:f][:fpkm] : '',
    ].join(',')
  end
  # Add header line
  lines.unshift("contig,ID male,locus male,FPKM male,ID female,locus female,FPKM female")
  # Write to output file
  file = File.open( PATH+"mergedM&FmauWholeGenome.csv","w" )
  file << lines.join("\n")
  file.close
end

main # Execute the code
