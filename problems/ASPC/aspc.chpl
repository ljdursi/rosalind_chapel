module ASPC {
  use BigInteger;
  param defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=1000000 : int;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n, m: int;
    var channel = stdin;
    if infile != defaultfilename {
      try! {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }
    }

    try! {
      var success = channel.read(n, m);
    }
    return (n, m);
  }

  proc n_choose_k_mod(n, k, modulus) : int {
    if k < 0 || k > n then
      return 0;

    if k == 0 || k == n then
      return 1;

    var num, den = new bigint(1);
    for i in 1..min(k, n-k) {
        num *= (n-i+1);
        den *= i;
    }
    return ((num / den) % modulus):int;
  }

  proc main() {
    try! {
      const (n, m) = params_from_file(infile, defaultfilename);
      var res = 0;
      for k in m..n do
        res = (res + n_choose_k_mod(n, k, modulus)) % modulus;
      writeln(res);
    }
  }
}
