module ORF {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const codonfile="codons.txt";

  // use stdin if filename == defaultfilename
  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.read(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line;
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  // read in the table of codons
  // because I'm too lazy to hardcode it
  proc readtables() throws {
    var allcodons : domain(string);
    var codons : [allcodons] string;
    var stop : domain(string);
    var start : domain(string);

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

    start += 'AUG';
    return (codons, start, stop);
  }

  proc reversecomplement(line: string) : string {
    var rc:string = "";
    const complements = ['A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'];

    for i in 1..line.length by -1 do
      rc += complements[line[i]];

    return rc;
  }

  proc transcribe(line: string) : string {
    return line.replace("T", "U");
  }

  proc translate(rnaseq, codons, start, stop) : domain(string) {
    var aminos = "";

    // generate the string of amino acids,
    // and identify start and stop locations
    var startlocs: [1..0] int;
    var stoplocs: [1..0] int;
    for i in 1..rnaseq.len-2 by 3 {
      var codon = rnaseq(i..i+3-1);
      if stop.member(codon) {
        stoplocs.push_back((i-1)/3+1);
        aminos += "@";
      } else {
        aminos += codons(codon);
        if start.member(codon) then
            startlocs.push_back((i-1)/3+1);
      }
    }
    
    var proteins : domain(string);
    var laststop = 0;
    for stop in stoplocs {
        for start in startlocs {
            if start > laststop && start <= stop then
                proteins.add(aminos[start..stop-1]) ;
        }
        laststop = stop;
    }
    return proteins;
  }

  proc allproteins(dnaseq, codons, start, stop) : domain(string) {
    var rnaseq = transcribe(dnaseq);
    var proteins: domain(string);
    for i in 1..3 do
        proteins = proteins + translate(rnaseq(i..), codons, start, stop);

    rnaseq = transcribe(reversecomplement(dnaseq));
    for i in 1..3 do
        proteins = proteins + translate(rnaseq(i..), codons, start, stop);

    return proteins;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    var tables = try! readtables();
    var codons = tables[1];
    var start = tables[2];
    var stop = tables[3];

    var sequences : [1..0] string;
    try! {
        for (seqlabel, sequence) in fasta_iterator(inchannel) {
          sequences.push_back(sequence);
        }

        for sequence in sequences {
            var proteins = allproteins(sequence, codons, start, stop);
            for protein in proteins do
                writeln(protein);
        }
    }
  }
}
