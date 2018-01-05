module SUBS {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var sequence: string = "";
    var query: string = "";

    try! {
      var channel = stdin;
      if infile != defaultfilename {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }

      var success = channel.read(sequence);
      success = channel.read(query);
    }

    return (sequence, query);
  }

  proc main() {
    var seq, query : string;
    try! {
      (seq, query) = params_from_file(infile, defaultfilename);
      const nseq = seq.len;
      const nquery = query.len;
        
      var locs : [1..0] int;
      for i in 1..nseq-nquery+1 do
          if seq(i..#nquery) == query then
            locs.push_back(i);
    
      writeln(locs);
    }
  }
}
