module EUBD {
  // Enumerating Unrooted Binary Trees
  // http://rosalind.info/problems/EUBT/
  param defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var line: string;
    var lines: [1..0] string;
    while (channel.readline(line)) do
      lines.push_back(line.strip());
  
    return lines;
  }

  record Edge {
    var node1idx : int;
    var node2idx : int;
  }

  record Node {
    var lab: string;
    var edgeidxs: [1..3] int;
    var nedges: int;
    var idx: int;

    proc add_edge(edgeidx) {
      nedges += 1;
      edgeidxs[nedges] = edgeidx;
    }

    proc remove_edge(edgeidx) {
      const (found, ei) = edgeidxs.find(edgeidx);

      if !found then
        return;

      edgeidxs[ei] = edgeidxs[nedges];
      nedges -= 1;
    }
  }

  class Tree {
    var nodes : [1..0] Node;
    var edges : [1..0] Edge;

    proc init(nodelabels: [] string, edges: [] 2*int) {
      this.complete();

      for lab in nodelabels do
        add_node(lab);

      for edge in edges do
        add_edge(edge[1], edge[2]);
    }

    proc init(other_tree: Tree) {
      this.complete();

      for edge in other_tree.edges do
        edges.push_back(edge);

      for node in other_tree.nodes do
        nodes.push_back(node);
    }

    proc add_node(nodelabel) {
      nodes.push_back(new Node(nodelabel, [-1, -1, -1], 0, 0));
      const nodeidx = nodes.size;
      nodes[nodeidx].idx = nodeidx;
      return nodeidx;
    }

    iter node_neighbours(nodeidx) {
      const node = nodes[nodeidx];
      for i in 1..node.nedges {
        const edgeidx = node.edgeidxs[i];
        const (n1, n2) = (edges[edgeidx].node1idx, edges[edgeidx].node2idx);
        const neighbour = if nodeidx == n1 then n2 else n1;
        yield neighbour;
      }
    }

    proc add_edge(node1, node2) {
      const edgenum = edges.size + 1;

      edges.push_back(new Edge(node1, node2));
      nodes[node1].add_edge(edgenum);
      nodes[node2].add_edge(edgenum);

      return edgenum;
    }

    proc split_edge(edgeidx, nodelabel) {
      const newnode = add_node(nodelabel);

      var (n1, n2) = (edges[edgeidx].node1idx, edges[edgeidx].node2idx);
      nodes[n2].remove_edge(edgeidx);
      edges[edgeidx].node2idx = newnode;
      nodes[newnode].add_edge(edgeidx);

      const newedge = add_edge(newnode, n2);
      
      return newnode;
    }

    proc insert_node_on_edge(edgeidx, nodelabel) {
      const internalnode = split_edge(edgeidx, "");
      const newnode = add_node(nodelabel);
      const newedge = add_edge(internalnode, newnode);
    }

    proc newick_dfs(nodeidx, visited) : string {
      var subtreestr = "";

      if visited[nodeidx] || nodeidx < 0 || nodeidx > nodes.size then
        return subtreestr;
  
      var open_neighbours : [1..0] int;
      var neigh: int;
      for neigh in this.node_neighbours(nodeidx) do
        if ! visited[neigh] then
          open_neighbours.push_back[neigh];
      const nneigh = open_neighbours.size;

      visited[nodeidx] = true;
      if open_neighbours.size == 0 then
        return nodes[nodeidx].lab;

      subtreestr = "(";
      for neigh in open_neighbours[1..nneigh-1] do
        subtreestr = subtreestr + newick_dfs(neigh, visited) + ",";

      subtreestr = subtreestr + newick_dfs(open_neighbours[nneigh], visited) + ")";
      subtreestr = subtreestr + nodes[nodeidx].lab;

      return subtreestr;
    }

    proc newick() {
      var visited : [1..nodes.size] bool = false;

      if nodes.size == 0 then
        return "";

      return newick_dfs(1, visited);
    }
  }

  proc enumerate_trees(const ref base_tree: Tree, newnodes: [?d] string) {
    if newnodes.size == 0 {
      writeln(base_tree.newick() + ";");
      return;
    }

    for edgeidx in 1..base_tree.edges.size {
      var subtree = new Tree(base_tree);
      subtree.insert_node_on_edge(edgeidx, newnodes[newnodes.domain.low]);
      enumerate_trees(subtree, newnodes[newnodes.domain.low+1..]);
    }
  }
  
  proc main() {
    try! {
      for line in lines_from_file(infile, defaultfilename) {
        const allnames = line.split();
        const tree = new Tree(allnames[1..2], [(1,2)]);
        enumerate_trees(tree, allnames[3..]);
      }
    }
  }
}
