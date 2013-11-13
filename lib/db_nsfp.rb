require 'rbbt-util'
require 'rbbt/workflow'

Workflow.require_workflow "Genomics"

module DbNSFP
  extend Workflow
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

  DbNSFP.claim DbNSFP.scores, :proc do |filename|
    organism = "Hsa/jan2013"

    require 'rbbt/sources/organism'

    files = DbNSFP.data.produce.glob('*variant*')

    transcript2protein = Organism.transcripts(organism).tsv :fields => ["Ensembl Protein ID"], :type => :single, :persist => true, :unnamed => true

    save_header = true
    Persist.persist_tsv("dbNSFP", nil, {}, {:file => filename, :persist => true, :prefix => "dbNSFP", :engine => TokyoCabinet::BDB, :serializer => :float_array}) do |database|
      files.each do |file|
        Log.info "Opening #{ file }"

        all_fields = TSV.parse_header(file).all_fields
        scores = all_fields.select{|f| f =~ /score/}

        if save_header
          database.fields = scores
          database.key_field = "Mutated Isoform"
          print_header = false
        end

        mutation_fields = %w(aaref aapos aaalt).collect{|f| all_fields.index f}
        transcript_field = all_fields.index "Ensembl_transcriptid"
        score_fields = scores.collect{|f| all_fields.index f}

        Open.read(file) do |line|
          next if line[0] == "#"
          parts = line.strip.split("\t").collect{|s| s == '.' ? "" : s }
          transcripts = parts[transcript_field].split ";"
          if transcripts.length == 1
            transcript = transcripts.first
            protein = transcript2protein[transcript]
            next if protein.nil? or protein.empty?

            mutation = parts.values_at(*mutation_fields) * ""
            next if mutation.empty?

            scores = parts.values_at(*score_fields)

            isoform = protein + ":" << mutation
            database[isoform] = scores.collect{|s| (s.nil? || s.empty?) ? -999 : s.to_f }
          else
            proteins = transcript2protein.values_at *transcripts
            next if proteins.compact.empty?

            m = parts.values_at(*mutation_fields)
            next if m.first.empty?
            mutations_zip = Misc.zip_fields(m.collect{|s| 
              p = s.split ";" 
              case p.length
              when 1
                [p.first] * proteins.length
              when proteins.length
                p
              else
                puts Hash[*all_fields.zip(parts).flatten].to_yaml
                raise "Number mismatch: #{proteins.length} proteins but #{p.length} parts"
              end
            })

            s = parts.values_at(*score_fields)
            scores_zip= Misc.zip_fields(s.collect{|s| 
              p = s.split ";" 
              case p.length
              when 0
                [nil] * proteins.length
              when proteins.length
                p
              else
                [p.first] * proteins.length
              #else
              #  puts Hash[*all_fields.zip(parts).flatten].to_yaml
              #  raise "Number mismatch: #{proteins.length} proteins but #{p.length} parts"
              end
            })

            transcripts.zip(proteins, mutations_zip, scores_zip).each do |transcript, protein,mutation_parts,scores|
              next if mutation_parts.nil? or scores.nil?
              next "Missing protein for transcript: #{ transcript }" if protein.nil? or protein.empty?
              isoform = protein + ":" << (mutation_parts * "")
              database[isoform] = scores.collect{|s| (s.nil? || s.empty?) ? -999 : s.to_f }
            end
          end
        end
      end

      database.type = :list
      database.namespace = organism
      database.cast = :to_f
      database.close
      database
    end

    nil
  end

  def self.database
    @@database ||= Persist.open_tokyocabinet(DbNSFP.scores.find, false, nil, TokyoCabinet::BDB)
  end

end
