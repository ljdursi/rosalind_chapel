module PPER {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=1000000 : uint;

  // use stdin if filename == defaultfilename
  proc param_from_file(filename: string, defaultfilename: string) throws {
    var n, m: uint;

    var channel = stdin;
    if infile != defaultfilename {
      try! {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }
    }

    try! {
      var success = channel.read(n);
      success = channel.read(m);
    }
    return (n, m);
  }

  proc partial_factorial_modulus(n, k, modulus) {
    if n <= 0 then
      return 0;

    var result : uint = 1;
    for i in n-k+1..n {
      result *= i;
      result %= modulus;
    }
    return result;
  }

  proc main() {
    try! {
      var params = param_from_file(infile, defaultfilename);
      writeln(partial_factorial_modulus(params[1], params[2], modulus));
    }
  }
}
