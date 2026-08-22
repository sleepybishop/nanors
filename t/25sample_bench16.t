use strict;
use warnings;
use Test::More;
use t::Util;

sub run_bench {
    my ($K, $N, $T) = @_;
    local $ENV{BENCH_MB} = $ENV{BENCH_MB} || 1;
    my ($stderr, $stdout) = run_prog("./t/00util/bench16 $K $N $T");
    return $stdout;
}

subtest "sample benchmark 16" => sub {
    my $T = 1280;
    my @targets = ("5-2", "7-2", "10-3", "20-5", "50-5", "100-10", "200-20", "500-25", "500-50");
    if ($ENV{CI_EMULATION}) {
        @targets = ("5-2", "7-2", "10-3", "20-5");
        $T = 256;
    }
    foreach (@targets) {
      my ($K, $N) = split("-", $_);
      my $resp = run_bench($K, $N, $T);

      my ($enc_mbps) = $resp =~ /encoded.*throughput: ([0-9\.]+)MB/;
      my ($dec_mbps) = $resp =~ /decoded.*throughput: ([0-9\.]+)MB/;

      diag "RS16_$_ T: $T enc/dec MBps: $enc_mbps/$dec_mbps";
      my $min_enc = $ENV{MIN_ENC_MBPS} || 0;
      my $min_dec = $ENV{MIN_DEC_MBPS} || 0;
      ok $enc_mbps > $min_enc, "encode throughput exceeds $min_enc MB/s for RS16_$_";
      ok $dec_mbps > $min_dec, "decode throughput exceeds $min_dec MB/s for RS16_$_";
    }
};

done_testing();
