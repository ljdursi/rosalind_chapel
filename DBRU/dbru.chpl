module DBRU {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  iter lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var line: string;
    while (channel.readline(line)) do
      yield line.strip();
  }

  proc reversecomplement(line: string) : string {
    var rc:string = "";
    const complements = ['A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'];

    for i in 1..line.length by -1 do
      rc += complements[line[i]];

    return rc;
  }

  proc main() {
    try! {
      var read: string;
      var readset: domain(string);

      for read in lines_from_file(infile, defaultfilename) {
        readset.add(read);        
        readset.add(reversecomplement(read));
      }
    
      for read in readset {
        const l = read.length;
        writeln("(", read(1..l-1), ", ", read(2..l), ")");
      }
    }
  }
}
