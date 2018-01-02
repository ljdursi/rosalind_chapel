module SEXL {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  iter reals_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var num: real;
    
    while (channel.read(num)) do
      yield num;
  }

  proc main() {
    try! {
      var cont = false;
      for q in reals_from_file(infile, defaultfilename) {
        if cont then
          write(" ");
        else
          cont = true;

        var p = 1.0-q;
        write(2.0*p*q);
      }
      writeln();
    }
  }
}
