module PPER {
  param defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=1000000 : uint;
  config const n=2 : uint;

  // use stdin if filename == defaultfilename
  proc param_from_file(filename: string, defaultfilename: string) throws {
    var k: uint;

    var channel = stdin;
    if infile != defaultfilename {
      try! {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }
    }

    try! {
      var success = channel.read(k);
    }
    return k;
  }

  proc n_to_the_k_mod(n, k, modulus) {
    if k <= 0 then
      return 1;

    var result : uint = 1;
    for i in 1..k {
      result *= n;
      result %= modulus;
    }
    return result;
  }

  proc main() {
    try! {
      const k = param_from_file(infile, defaultfilename);
      writeln(n_to_the_k_mod(n, k, modulus));
    }
  }
}
