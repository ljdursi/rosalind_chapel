module CUNR {
  const defaultfilename = "none";
  config const infile = defaultfilename;
  config const modulo = 1000000;

  // use stdin if filename == defaultfilename
  iter ints_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var num: int;
    
    while (channel.readln(num)) do
      yield num;
  }

  proc double_factorial_mod(n, mod) {
    var result = 1;
    for i in 1..n by -2 do
      result = (result * i) % mod;
    
    return result;
  }

  proc main() {
    try! {
      for n in ints_from_file(infile, defaultfilename) do
        writeln(double_factorial_mod(2*n-5, modulo));
    }
  }
}
