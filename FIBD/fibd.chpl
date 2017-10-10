module FIB {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;
    var m: int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n, m);
    return (n, m);
  }

  proc mortal_fib(ngen:int, lifespan:int): uint {
    var pop_structure: [1..lifespan] uint;
    pop_structure = 0;
    pop_structure[1] = 1;

    for i in 2..ngen {
      var offspring = + reduce pop_structure[2..];
      pop_structure[2..] = pop_structure[1..lifespan-1];
      pop_structure[1] = offspring;
    }

    return + reduce pop_structure;
  }

  proc main() {
    var n:int, m:int;
    try! {
      (n, m) = params_from_file(infile, defaultfilename);
      var result = mortal_fib(n, m);

      writeln(result);
    }
  }
}
