module DNA {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // read file into a single string 
  // use stdin if filename == defaultfilename
  proc string_from_file(filename: string, defaultfilename: string) : string {
    var text: string = "";
    var line: string = "";

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.read(line)) do 
      text = text + line;

    return text;
  }

  proc count_chars(line: string) {
    var chars : domain(string);
    var counts: [chars] int;
    for c in line {
      counts[c] += 1;
    }
    return counts;
  }

  proc main() {
    var sequence = string_from_file(infile, defaultfilename);
    var counts = count_chars(sequence);

    writeln(counts['A'], " ", counts['C'], " ", counts['G'], " ", counts['T']);
  }
}
