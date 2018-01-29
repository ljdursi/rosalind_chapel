module FOUN {
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;
    var probs: [1..0] real;
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.readln(n);
    channel.readline(line);
    line.strip();
    for item in line.split() do
      probs.push_back(item:real);

    return (n, probs);
  }

  proc main() {
    try! {
      const (npop, probs) = params_from_file(infile, defaultfilename);
      writeln(npop * probs);
    }
  }
}
