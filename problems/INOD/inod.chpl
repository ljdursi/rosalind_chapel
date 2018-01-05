module INOD {
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n);
    return n;
  }

  proc main() {
    try! {
      var nleaves = params_from_file(infile, defaultfilename);
      writeln(nleaves - 2);
    }
  }
}
