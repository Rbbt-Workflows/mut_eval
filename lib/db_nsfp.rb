require 'rbbt-util'
$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))

module DbNSFP
  extend Resource
  self.subdir = 'var/dbNSFP'

  DbNSFP.claim DbNSFP.data, :proc do |directory|
    url = "http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.1.zip"
    Misc.in_dir(directory) do
      FileUtils.mkdir_p '.source'
      `wget '#{url}' -c -O .source/pkg.zip && cd .source && unzip pkg.zip && find . -name '*variant*' | xargs -I '{}' mv '{}' ../ && cd .. && rm -Rf .source *.zip`
    end
    nil
  end

  def self.database
    @@database ||= begin
                     sharder_function = Proc.new do |key|
                       key[(13..14)]
                     end

                     Persist.persist_tsv("dbNSFP", nil, {}, 
                                         :file => DbNSFP.scores_shard.find, :persist => true, :prefix => "dbNSFP", :serializer => :float_array,
                                         :shard_function => sharder_function) do |sharder|

                       require 'rbbt/sources/organism'

                       organism = "Hsa/jan2013"

                       files = DbNSFP.data.produce.glob('*variant*')

                       transcript2protein = Organism.transcripts(organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true, :unnamed => true

                       save_header = true
                       TSV.traverse files.sort, :bar => "DbNSFP files" do |file|
                         all_fields = TSV.parse_header(file).all_fields
                         scores = all_fields.select{|f| f =~ /score/}

                         if save_header
                           sharder.fields = scores
                           sharder.key_field = "Mutated Isoform"
                           print_header = false
                         end

                         mutation_fields = %w(aaref aapos aaalt).collect{|f| all_fields.index f}
                         transcript_field = all_fields.index "Ensembl_transcriptid"
                         score_fields = scores.collect{|f| all_fields.index f}

                         sharder.write
                         TSV.traverse file, :type => :array, :bar => File.basename(file), :into => sharder do |line|
                           next if line[0] == "#"
                           parts = line.strip.split("\t",-1)
                           transcripts = parts[transcript_field].split ";"

                           if transcripts.length == 1
                             transcript = transcripts.first
                             protein = transcript2protein[transcript]
                             next if protein.nil? or protein.empty?

                             mutation_parts = parts.values_at(*mutation_fields)
                             next if mutation_parts[1] == "-1"

                             scores = parts.values_at(*score_fields)

                             isoform = protein + ":" << mutation_parts * ""
                             values = scores.collect{|s| s == '.' ? -999 : s.to_f }

                             [isoform, values]
                           else
                             proteins = transcript2protein.values_at *transcripts
                             next if proteins.compact.empty?

                             mutation_parts = parts.values_at(*mutation_fields)
                             next if mutation_parts[1] == "-1"

                             if mutation_parts[1].index ";"
                               mutations_zip = mutation_parts[1].split(";").collect{|pos| [m[0],pos,m[2]] }
                             else
                               mutations_zip = [mutation_parts] * proteins.length
                             end

                             s = parts.values_at(*score_fields)

                             scores_zip = [s] * proteins.length

                             results = [].extend TSV::MultipleResult
                             transcripts.each_with_index do |transcript,i|
                               protein = proteins[i]
                               next if protein.nil? or protein.empty?
                               isoform = protein + ":" << (mutations_zip[i] * "")
                               values = scores_zip[i].collect{|s| s == '.' ? -999 : s.to_f }
                               results << [isoform, values]
                             end

                             results
                           end
                         end
                       end # traverse files

                       sharder.read; 
                     end # persist
                   end # end
  end

end
