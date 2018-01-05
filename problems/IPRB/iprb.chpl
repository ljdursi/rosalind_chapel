module IPRB {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var k, m, n: uint;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(k, m, n);
    return (k, m, n);
  }

  proc generalized_fib(ngen:int, off_per_gen:int): int {
    var last = 0, current = 1;
    for i in 2..ngen do
      (last, current) = (current, current + off_per_gen*last );

    return current;
  }

  proc main() {
    var k, m, n : uint;
    try! {
      (k, m, n) = params_from_file(infile, defaultfilename);
      var tot = k + m + n;
      var p = (1.0/tot)*(k + m/(tot-1.0)*(k + 0.75*(m-1.0) + 0.5*n) + n/(tot-1.0)*(k + 0.5*m));

      writeln(p);
    }
  }
}
