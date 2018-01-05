module EVAL {
  use Math;
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

  proc string_log_prob(sequence, gcprob) {
    var nbases = ["A" => 0, "C" => 0, "G" => 0, "T" => 0];
    for base in sequence do
        nbases[base] += 1;

    const loggprob = log(gcprob/2.0);
    const logaprob = log((1-gcprob)/2.0);

    var logprob = nbases["G"]*loggprob + nbases["C"]*loggprob;
    logprob += nbases["A"]*logaprob + nbases["T"]*logaprob;

    return logprob;
  }

  proc expectation(sequence, nchances, gcprob) {
    const log_nexpect = log(nchances) + string_log_prob(sequence, gcprob);
    return exp(log_nexpect);
  }

  proc main() {
    try! {
      var channel = input_channel(infile, defaultfilename);
      var seq, line : string;

      channel.readline(line);
      line.strip();
      var n = line : int;

      channel.readline(line);
      seq = line.strip();

      const nchances = n - seq.length + 1;

      channel.readline(line);
      line = line.strip();
      var nexpect : [1..0] real;

      for item in line.split() do
        nexpect.push_back(expectation(seq, nchances, item:real));

      writeln(nexpect);
    }
  }
}
