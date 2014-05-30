$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'db_nsfp'

class TestClass < Test::Unit::TestCase
  def _test_kch
    database = Persist.open_tokyocabinet(DbNSFP.scores_kch.produce.find, false, nil, "kch")
  end

  def test_shard
   DbNSFP.scores_shard.produce
  end
end

