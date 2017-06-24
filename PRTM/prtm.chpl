module RNA {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const weightfile="weights.txt";

  proc readtable() {
    var aminos : domain(string);
    var weights : [aminos] real;

    var file = open(weightfile, iomode.r);
    var channel = file.reader();

    var amino : string;
    var weight : real;
    while (channel.readln(amino, weight)) do
      weights[amino] = weight;

    return weights;
  }

  // use stdin if filename == defaultfilename
  proc string_from_file(filename: string, defaultfilename: string) : string {
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

  proc totalweight(protein : string, weights) : real {
    var weight = 0.0;
    for amino in protein {
      weight += weights[amino]; 
    }
    return weight;
  }

  proc main() {
    var weights = readtable();
    var protein = string_from_file(infile, defaultfilename);
    var result = totalweight(protein, weights);
     
    writef("%dr\n", result);
  }
}
