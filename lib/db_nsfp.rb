require 'rbbt-util'
$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))

module DbNSFP
  extend Resource
  self.subdir = 'var/dbNSFP'

  MI_SHARD_FUNCTION = Proc.new do |key|
    key[(13..14)]
  end

  GM_SHARD_FUNCTION = Proc.new do |key|
    key[0..key.index(":")-1]
  end

  MI_SCORES = %w(SIFT_score SIFT_converted_rankscore SIFT_pred Polyphen2_HDIV_score
                 Polyphen2_HDIV_rankscore Polyphen2_HDIV_pred Polyphen2_HVAR_score
                 Polyphen2_HVAR_rankscore Polyphen2_HVAR_pred LRT_score LRT_converted_rankscore
                 LRT_pred MutationTaster_score MutationTaster_converted_rankscore
                 MutationTaster_pred MutationAssessor_score MutationAssessor_rankscore
                 MutationAssessor_pred FATHMM_score FATHMM_rankscore FATHMM_pred RadialSVM_score
                 RadialSVM_rankscore RadialSVM_pred LR_score LR_rankscore LR_pred
                 Reliability_index VEST3_score VEST3_rankscore)

  def self.organism
    "Hsa/dec2013"
  end

  DbNSFP.claim DbNSFP.data, :proc do |directory|
    url = "http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.5.zip"
    Misc.in_dir(directory) do
      FileUtils.mkdir_p '.source'
      `wget '#{url}' -c -O .source/pkg.zip && cd .source && unzip pkg.zip && find . -name '*variant*' | xargs -I '{}' mv '{}' ../ && cd .. && rm -Rf .source *.zip`
    end
    nil
  end

  DbNSFP.claim DbNSFP.mutations, :proc do |filename|
    Misc.sensiblewrite(filename) do |f|
      DbNSFP.data.glob('*variant*').each do |file|
        next unless file =~ /chr19/
        TSV.traverse file, :type => :array, :bar => file do |line|
          line = line.strip.gsub(/\t\.\t/, "\t\t")
          if line =~ /^#/
            f.puts("#: :type=:list#:namespace=#{organism}")
            f.puts("#Genomic Mutation" << "\t" << line[1..-1])
          else
            chr, pos, ref, alt = line.split("\t")
            f.puts(([chr, pos, alt] * ":" ) << "\t" << line)
          end
        end
      end
    end
  end

  DbNSFP.claim DbNSFP.data.readme, :url, "http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.5.readme.txt"

  def self.database
    @@database||= begin
                     Persist.persist_tsv("dbNSFP", nil, {}, :persist => true, :update => false,
                                         :file => DbNSFP.scores_shard.find, :prefix => "dbNSFP", :serializer => :float_array, :engine => "BDB",
                                         :shard_function => MI_SHARD_FUNCTION) do |sharder|

                       require 'rbbt/sources/organism'

                       organism = "Hsa/jan2013"

                       files = DbNSFP.data.produce.glob('*variant*')#.select{|f| f =~ /22/}

                       transcript2protein = Organism.transcripts(organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true, :unnamed => true

                       save_header = true
                       TSV.traverse files.sort, :bar => "DbNSFP files" do |file|
                         all_fields = TSV.parse_header(file).all_fields
                         scores = all_fields[23..-1]

                         if save_header
                           sharder.fields = scores
                           sharder.key_field = "Mutated Isoform"
                           print_header = false
                         end

                         mutation_fields = %w(aaref aapos aaalt).collect{|f| all_fields.index f}
                         transcript_field = all_fields.index "Ensembl_transcriptid"
                         score_fields = scores.collect{|f| all_fields.index f}

                         TSV.traverse file, :type => :array, :bar => File.basename(file) do |line|
                           next if line[0] == "#"

                           parts = line.strip.split("\t",-1)
                           transcripts = parts[transcript_field].split ";"

                           res = if transcripts.length == 1
                             transcript = transcripts.first
                             protein = transcript2protein[transcript]
                             next if protein.nil? or protein.empty?

                             mutation_parts = parts.values_at(*mutation_fields)
                             next if mutation_parts[1] == "-1"

                             scores = parts.values_at(*score_fields)

                             isoform = protein + ":" << mutation_parts * ""
                             values = scores.collect{|s| s == '.' ? -999 : s.to_f }

                             sharder[isoform] = values
                             #[isoform, values]
                           else
                             proteins = transcript2protein.values_at *transcripts
                             next if proteins.compact.empty?

                             mutation_parts = parts.values_at(*mutation_fields)
                             next if mutation_parts[1] == "-1"

                             if mutation_parts[1].index ";"
                               mutations_zip = mutation_parts[1].split(";").collect{|pos| [mutation_parts[0],pos,mutation_parts[2]] }
                             else
                               mutations_zip = [mutation_parts] * proteins.length
                             end

                             s = parts.values_at(*score_fields)

                             scores_zip = [s] * proteins.length

                             #results = [].extend MultipleResult

                             transcripts.each_with_index do |transcript,i|
                               protein = proteins[i]
                               next if protein.nil? or protein.empty?
                               isoform = protein + ":" << (mutations_zip[i] * "")
                               values = scores_zip[i].collect{|s| s == '.' ? -999 : s.to_f }
                               sharder[isoform] = values
                               #results << [isoform, values]
                             end

                             #results
                           end
                         end
                       end # traverse files

                     end # persist
                   end # end
  end
end
