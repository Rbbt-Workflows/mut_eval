$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'sharder'

class TestClass < Test::Unit::TestCase
  def test_shard
    TmpFile.with_file do |dir|
      sharder = Sharder.new dir, true, :float_array, 'kch' do |key|
        key[-1]
      end

      1000.times do 
        v = rand(100000)
        sharder[v.to_s] = [v, v*2]
      end
      sharder.read

      ddd sharder.keys
    end
  end
end

