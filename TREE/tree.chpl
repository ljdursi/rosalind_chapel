module TREE {
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;
    var starts: [1..0] int;
    var ends: [1..0] int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n);
    while success {
        var start, end : int;
        success = channel.read(start, end);
        if success {
            starts.push_back(start);
            ends.push_back(end);
        } 
    }
    return (n, starts, ends);
  }

  proc main() {
    try! {
      var (nnodes, edgestarts, edgeends) = params_from_file(infile, defaultfilename);
      writeln(nnodes - 1 - edgestarts.size);
    }
  }
}
