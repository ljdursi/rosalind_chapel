module FIB {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;
    var k: int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n, k);
    return (n, k);
  }

  proc generalized_fib(ngen:int, off_per_gen:int): int {
    var last = 0, current = 1;
    for i in 2..ngen do
      (last, current) = (current, current + off_per_gen*last );

    return current;
  }

  proc main() {
    var n:int, k:int;
    try! {
      (n, k) = params_from_file(infile, defaultfilename);
      var result = generalized_fib(n, k);

      writeln(result);
    }
  }
}
