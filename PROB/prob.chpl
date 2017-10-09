module PROB {
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  proc stringprobs(sequence, gcprobs) {
    var nbases = ["A" => 0, "C" => 0, "G" => 0, "T" => 0];
    for base in sequence do
        nbases[base] += 1;

    var logprobs : [1..0] real;
    for gcprob in gcprobs {
        const loggprob = log10(gcprob/2.0);
        const logaprob = log10((1-gcprob)/2.0);

        var logprob = nbases["G"]*loggprob + nbases["C"]*loggprob;
        logprob += nbases["A"]*logaprob + nbases["T"]*logaprob;

        logprobs.push_back(logprob);
    }
    return logprobs;
  }

  proc main() {
    try! {
      var channel = input_channel(infile, defaultfilename);
      var seq, values : string;

      channel.readline(seq);
      channel.readline(values);
      var gcprobs : [1..0] real;
      for value in values.split() do
        gcprobs.push_back(value : real);

      writeln(stringprobs(seq, gcprobs));
    }
  }
}
