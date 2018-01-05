module PROT {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const codonfile="codons.txt";

  proc readtables() throws {
    var allcodons : domain(string);
    var codons : [allcodons] string;
    var stop : domain(string);

    try! {
      var file = open(codonfile, iomode.r);
      var channel = file.reader();

      var rna : string;
      var amino : string;
      while (channel.readln(rna, amino)) {
        if amino == "Stop" then
          stop += rna;
        else
          codons[rna] = amino;
      }
    }

    return (codons, stop);
  }

  // use stdin if filename == defaultfilename
  proc string_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";
    var line: string = "";

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readln(line)) do 
      text += line;

    return text;
  }

  proc transcribe(rnaseq, codons, stop) : string {
    var protein = "";

    for i in 1..rnaseq.len-2 by 3 {
      var codon = rnaseq(i..i+3-1);
      if stop.member(codon) then
        return protein;
      else
        protein += codons[codon];
    }
    return protein;
  }

  proc main() {
    try! {
      var tables = readtables();
      var codons = tables[1];
      var stop = tables[2];

      var rna = string_from_file(infile, defaultfilename);
      var result = transcribe(rna, codons, stop);
     
      writeln(result);
    }
  }
}
