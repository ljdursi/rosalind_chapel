module SPLC {
  use Regexp;
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
        sequence += line.strip();
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

  proc transcribe(line: string) : string {
    return line.replace("T", "U");
  }

  proc translate(rnaseq, codons, start, stop) {
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
    
    var proteins : [1..0] string;
    var laststop = 0;
    for stop in stoplocs {
        for start in startlocs {
            if start > laststop && start <= stop then
                proteins.push_back(aminos[start..stop-1]) ;
        }
        laststop = stop;
    }
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

        var dna = sequences[1];
        var introns = sequences[2..];

        for intron in introns {
            var intronre = compile(intron);
            dna = intronre.sub("", dna);
        }
        const rna = transcribe(dna);
        const proteins = translate(rna, codons, start, stop);
        writeln(proteins[1]);
    }
  }
}
