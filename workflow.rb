require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/persist/tsv'
#require 'rbbt/mutation/mutation_assessor'
#require 'rbbt/mutation/polyphen'
#require 'rbbt/mutation/transFIC'
#require 'rbbt/mutation/sift'
#require 'rbbt/mutation/snps_and_go'

Workflow.require_workflow "Genomics"

require 'rbbt/entity/protein'
require 'rbbt/util/R'

#require 'db_nsfp'

module MutEval
  extend Workflow

  CACHE_DIR = Rbbt.var.cache.MutEval.find

  CACHES = {
    :missing_predictions => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'missing_predictions'), false, :string),
    :mutation_assessor => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'mutation_assessor'), false, :string),
    :sift => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'sift'), false, :string),
    :polyphen => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'polyphen'), false, :string),
    :snps_and_go => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'snps_and_go'), false, :string),
    :transFIC => Persist.open_tokyocabinet(File.join(CACHE_DIR, 'transFIC'), false, :string),
  }

  CACHES.values.each{|db| db.close}

  helper :get_cache do |method, mutations|
    tsv = TSV.setup({})

    cache = CACHES[method]
    missing = []

    cache.read_and_close do
      mutations.each do |mutation|
        if cache.include? mutation
          tsv[mutation] = cache[mutation].split("\t") 
        else
          missing << mutation
        end
      end
    end

    Log.debug("Cache for #{ method } found #{tsv.length} out of #{mutations.length}; missing #{missing.length}")
    [tsv, missing]
  end

  helper :add_cache do |method, tsv|
    cache = CACHES[method]
    cache.write_and_close do
      tsv.each do |mutation, values|
        cache[mutation] = values * "\t"
      end
    end
  end

  helper :get_missing do |method, mutations|
    tsv = TSV.setup({})

    cache = CACHES[:missing_predictions]
    missing = []

    cache.read_and_close do
      mutations.each do |mutation|
        if cache.include? mutation and cache[mutation].include? method.to_s
          missing << mutation
        end
      end
    end

    Log.debug("Mutations known to have missing predictions by #{ method }: #{missing.length} out of #{mutations.length}")

    missing
  end


  helper :add_missing do |method, mutations, tsv|
    cache = CACHES[:missing_predictions]
    cache.write_and_close do
      (mutations - tsv.keys).each do |mutation|
        if cache.include? mutation
          cache[mutation] << method.to_s
        else
          cache[mutation] = [method.to_s]
        end
      end
    end
  end

  helper :get_dbNSFP do |mutations,method|
    mutations = mutations.select{|m| m =~ /:([A-Z])\d+([A-Z])/ and $1 != $2 }.uniq.sort
    field = case method.to_s
            when "all"
              nil
            when "mutation_assessor"
              "MutationAssessor_score_converted"
            when "sift"
              "SIFT_score_converted"
            when "polyphen"
              "Polyphen2_HDIV_score"
            when "ltr"
              "LRT_score_converted"
            when "mutation_tasser"
              "MutationTaster_score_converted"
            when "fathmm"
              "FATHMM_score_converted"
            end
    Log.low "Querying dbNSFP (#{method || "all"}) with #{mutations.length} mutations"
    database = DbNSFP.database
    database.read_and_close do
      if field.nil?
        database.select(:key => mutations)
      else
        database.select(:key => mutations).slice(field)
      end
    end
  end

  input :mutations, :array, "Mutated Isoforms", nil
  task :dbNSFP => :tsv do |mutations|
    database = DbNSFP.database
    database.unnamed = true
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => database.fields, :type => :list, :cast => :to_f
    dumper.init
    TSV.traverse mutations, :into => dumper, :bar => true, :type => :array do |mutation|
      p = database[mutation]
      next if p.nil?
      p.collect!{|v| v == -999 ? nil : v }
      [mutation, p]
    end
  end
  export_synchronous :dbNSFP
end
