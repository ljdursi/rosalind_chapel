module INDC {
  use BigInteger;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc param_from_file(filename: string, defaultfilename: string) throws {
    var n: uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n);
    return n;
  }

  proc n_choose_ks(n: uint) {
    // pascal's triangle
    var triangle : [0..n, 0..n] bigint;
    for i in 0..n do
        triangle[i, 0] = new bigint(1);

    for i in 1..n do
      for j in 1..i do
        triangle[i, j] = triangle[(i-1), (j-1)] + triangle[(i-1), (j)];

    return triangle[n, 0..n-1] : [1..n] bigint;
  }

  proc binom_sfd(n: uint) {
    var ncks = n_choose_ks(n);
    var prefix = -1.0*n*log10(2.0);
    var sum = new bigint(0);
    var tot = new bigint(1);
    tot <<= n;

    var sf : [1..n] real;
    for i in 1..n do {
        sum += ncks[i];
        if i < n/2 then
          sf[i] = log1p(-sum:real/tot:real)/log(10.0);
        else
          sf[i] = log10((tot-sum):real) + prefix;
    }
    return sf;
  }

  proc main() {
    var k:uint, n:uint;
    try! {
      var n = param_from_file(infile, defaultfilename);
      writeln(binom_sfd(2*n));
    }
  }
}
