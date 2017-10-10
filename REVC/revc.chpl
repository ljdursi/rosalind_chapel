module REVC {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // read file into a single string 
  // use stdin if filename == defaultfilename
  proc string_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";

    try! {
      var line: string = "";
      var channel = stdin;
      if infile != defaultfilename {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }

      while (channel.read(line)) do 
        text = text + line;
    }

    return text;
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
      var dna = string_from_file(infile, defaultfilename);
      var rc = reversecomplement(dna);

      writeln(rc);
    }
  }
}
