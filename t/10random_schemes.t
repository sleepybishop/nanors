use strict;
use warnings;
use Test::More;
use t::Util;
use File::Temp qw(tempfile tempdir);
use Digest::MD5 qw(md5_hex);
use Data::Dumper;

sub run_test {
    return run_prog("./t/00util/test");
}

subtest "random schemes" => sub {
    my ($idx, $iters) = (0, 200);
    while ($idx < $iters) {
        my $resp = run_test($_);
        like $resp, "/===OK===/", "iter $idx";
        $idx++;
    }
};
        
done_testing();
