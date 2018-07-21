module MULT {
  // Multiple Alignment
  // http://rosalind.info/problems/mult/
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

  iter fasta_iterator(inchannel) throws {
    // iterate over fasta sequences
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.readline(line)) {
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
    if started then
      yield (seqlabel, sequence);
  }

  proc delta_score(align_bases) {
    // Given a column of the alignment, eg ['A','A','T','-'],
    // return the change in score - -1*number of pairwise mismatches
    var score = 0;
    for i in 1..align_bases.size do
      for j in i+1..align_bases.size do
        if align_bases[i] != align_bases[j] then
          score -= 1;

    return score;
  }

  proc delta_idxs_from_move(move, param dimension) {
    // given a move, calculate the shift in indexes the move was from
    var delta_idxs : dimension*int;

    var partial_moves : [1..dimension] int;
    for i in 1..dimension do 
      partial_moves[i] = (2**(i-1));

    var pmove = -move;
    for i in 1..dimension by -1 do
        if pmove / partial_moves[i] == 1 {
          pmove -= partial_moves[i];
          delta_idxs[i] = -1;
        }

    return delta_idxs;
  }

  proc move_from_delta_idxs(delta_idxs) {
    // given a shift in indexes, calculate an integer move
    var scale = 1;
    var move = 0;
    for delta in delta_idxs {
      move += delta*scale;
      scale *= 2;
    }
    return move;
  }

  proc align_bases(seqs, delta_idxs, idxs) {
    // given a the sequences, the current position, and
    // the shift in indexes from this move, generate the
    // column in the alignment
    const dimension = seqs.size;
    var bases : [1..dimension] string;
    for d in 1..dimension {
      if delta_idxs[d] == 0 then
        bases[d] = '-';
      else {
        bases[d] = seqs[d](idxs[d]);
      }
    }
    return bases;
  }

  proc multiple_alignment(sequences) {
    param rank = 4;  // this needs to be hardcoded
    var mtx_sizes: rank*range;
    var stencil_sizes: rank*range;

    for d in 1..rank {
      mtx_sizes(d) = 0..sequences[d].length;
      stencil_sizes(d) = -1..0;
    }

    var mtx_domain: domain(rank) = mtx_sizes;
    var stencil: domain(rank) = stencil_sizes;

    var score, moves : [mtx_domain] int;

    for idx in mtx_domain {
      if + reduce idx == 0 then
        continue;
      score[idx] = -999999;

      for delta in stencil {
        const prev_idx = idx + delta;
        if min reduce prev_idx < 0 then
          continue;

        const bases = align_bases(sequences, delta, idx);
        const trial_score = score[prev_idx] + delta_score(bases);
        const trial_move = move_from_delta_idxs(delta);

        if trial_score > score[idx] {
          score[idx] = trial_score;
          moves[idx] = trial_move;
        }
      }
    }

    // backtrack
    var idx = mtx_domain.high;
    const dist = score[idx];
    var augmenteds : [1..rank] string;

    while + reduce idx > 0 {
      const delta = delta_idxs_from_move(moves[idx], rank);
      var prev_idx = idx + delta;
      const bases = align_bases(sequences, delta, idx);

      for d in 1..rank do 
        augmenteds[d] = bases[d] + augmenteds[d];
      
      idx = prev_idx;
    }

    return (dist, augmenteds);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var (d, augs) = multiple_alignment(seqs);
      writeln(d);
      for aug in augs do
        writeln(aug);
    }
  }
}
