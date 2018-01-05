module MRNA {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const codonfile="codons.txt";
  config const modulus = 1000000 : uint;

  proc readtables() throws {
    var aminoacids : domain(string);
    var codons : [aminoacids] uint;
    var nstops = 0 : uint;

    try! {
      var file = open(codonfile, iomode.r);
      var channel = file.reader();

      var rna : string;
      var amino : string;
      while (channel.readln(rna, amino)) {
        if amino == "Stop" then
          nstops += 1;
        else
          codons[amino] += 1;
      }

      return (codons, nstops);
    }
  }

  // use stdin if filename == defaultfilename
  proc string_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";
    var line: string = "";

    try! {
      var channel = stdin;
      if infile != defaultfilename {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }

      while (channel.readln(line)) do 
        text += line;

      return text;
    }
  }

  proc ncombinations(protein, ncodons, nstops) : uint {
    var ncomb = nstops : uint;
    for amino in protein {
      ncomb *= ncodons[amino];
      ncomb = ncomb % modulus;
    }
    return ncomb;
  }

  proc main() {
    try!{ 
      var (codons, nstops) = readtables();
      var protein = string_from_file(infile, defaultfilename);

      var result = ncombinations(protein, codons, nstops);
      writeln(result);
    }
  }
}
