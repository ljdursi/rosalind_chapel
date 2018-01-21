module PDPL {
  use Sort;
  const defaultfilename="none";
  config const infile=defaultfilename;

  iter ints_from_file(filename: string, defaultfilename: string) throws {
    var line: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    
    while (channel.readline(line)) {
      const sline = line.strip();
      const items = sline.split();
      var ints : [1..#items.size] int;
      for i in 1..#items.size do
        ints[i] = items[i]:int;

      sort(ints);
      yield ints;
    }
  }

  proc update_sorted_array(curvals, newval) {
    var newvals = curvals[1..curvals.size];
    newvals.push_back(newval);
    sort(newvals);
    return newvals;
  }

  proc update_diffs_with_new_value(diffs, values, newval) {
    var newdiffs = abs(values - newval);
    var updateddiffs = diffs[1..diffs.size];

    for newdiff in newdiffs {
      var (success, idx) = updateddiffs.find(newdiff);
      if !success then
        return (false, updateddiffs);
      else
        updateddiffs.remove(idx);
    }

    return (true, updateddiffs);
  }

  record solution {
    var success: bool;
    var D : domain(1);
    var values: [D] int;
  }

  proc solve(values, diffs) : solution {
    if (diffs.size == 0) then
      return new solution(true, values.domain, values);

    const nextdiff = diffs[diffs.size];
    const leftval = values[values.size] - nextdiff;
    const rightval = values[1] + nextdiff;

    var (lsuccess, leftdiffs) = update_diffs_with_new_value(diffs, values, leftval);
    if lsuccess {
      var leftsoln = solve(update_sorted_array(values, leftval), leftdiffs);
      if leftsoln.success then
        return leftsoln;
    }

    var (rsuccess, rightdiffs) = update_diffs_with_new_value(diffs, values, rightval);
    if !rsuccess then
      return new solution(false, values.domain, values);

    return solve(update_sorted_array(values, rightval), rightdiffs);
  }

  proc infer_set_from_multiset(dmultiset) {
    const ndiffs = dmultiset.size;
    var values = [0, dmultiset[ndiffs]];
    var diffs = dmultiset[1..ndiffs-1];

    return solve(values, diffs);
  }

  proc main() {
    try! {
      for diffmultiset in ints_from_file(infile, defaultfilename) {
        var soln = infer_set_from_multiset(diffmultiset);
        writeln(soln.values);
      }
    }
  }
}
