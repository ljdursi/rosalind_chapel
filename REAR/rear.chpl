module REAR {
  use List;
  const defaultfilename="none";
  config const infile=defaultfilename;

  proc strings_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";
    var line: string;
    var lines: [1..0] string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.readline(line)) {
      lines.push_back(line.strip());
    }

    return lines;
  }


  proc permutation_from_line(line) {
    const items = line.split();
    var permutation : [1..#items.size] int;
    for i in 1..#items.size do
        permutation[i] = items[i]:int;
    return permutation;
  }

  proc inverse_permutation(permutation) {
    const n = permutation.size;
    var inverse : [1..#n] int;
    for i in 1..n do
        inverse[permutation[i]] = i;
    return inverse;
  }

  proc apply_permutation(permutation, series) {
    const n = series.size;
    var result : [1..#n] int;
    for i in 1..#n do
        result[i] = permutation[series[i]]; 
    return result;
  }

  proc adjacent(a:int, b:int) : bool{
    return (a+1 == b) || (a == b+1);
  }

  proc breakpoints(series) {
    var bps : [1..0] int;
    const n = series.size;
    for i in 2..n do
      if !adjacent(series[i-1], series[i]) then
        bps.push_back(i);

    return bps;
  }

  proc potential_reversals(series, bkpts) {
    const nbkpts = bkpts.size;
    var max_new_adj = -1;
    var potentials : [1..0] (2*int);

    for i in 2..nbkpts {
      for j in 1..i-1 {
        if adjacent(bkpts[i], bkpts[j]) then
          continue;
        
        var new_l_adj = adjacent(series[bkpts[i]-1], series[bkpts[j]-1]);
        var new_r_adj = adjacent(series[bkpts[i]], series[bkpts[j]]);
        var nnew_adj = new_l_adj : int + new_r_adj : int;

        if nnew_adj > max_new_adj {
          potentials.clear();
          max_new_adj = nnew_adj;
        }
        if nnew_adj == max_new_adj then
          potentials.push_back((bkpts[j], bkpts[i])); 
      }
    }
    return potentials;
  }

  proc exhaustive_reversal_distance(sequence, best) : int {
    var bkpts = breakpoints(sequence);
    const n = bkpts.size;
    if n == 0 then
      return 0;

    // short circuit if there's already a better solution than this
    if best == 0 then
      return 2*n;

    var moves = potential_reversals(sequence, bkpts);
    var shortest_dist = 2*n;
    for (leftbp, rightbp) in moves {
      var newsequence = sequence;
      for i in leftbp..rightbp-1 do
        newsequence[i] = sequence[rightbp-1-(i-leftbp)];
      const dist = exhaustive_reversal_distance(newsequence, best-1);
      if dist < shortest_dist then
        shortest_dist = dist;
      if shortest_dist <= ((n+1)/2):int then
        break;
    }
    
    return 1 + shortest_dist;
  }

  proc main() {
    try! {
      var lines = strings_from_file(infile, defaultfilename);
      var distances : [1..0] int;
      for (line1, line2) in zip(lines[1.. by 3], lines[2.. by 3]) {
        var perm1 = permutation_from_line(line1);
        var perm2 = permutation_from_line(line2);
        const n = perm1.size;

        var normalized_perm = apply_permutation(inverse_permutation(perm1), perm2);
        var extended_perm : [1..n+2] int;
        extended_perm[1] = 0;
        extended_perm[2..#n] = normalized_perm;
        extended_perm[n+2] = n+1;
        distances.push_back(exhaustive_reversal_distance(extended_perm, n+1));
      }
      writeln(distances);
    }
  }
}
