sfdp-lite
=========

sfdp-lite is a simple graph layout tool that relies on the Graphviz sfdp algorithm to compute the layout of large graphs in Matrix Market format.

<img src="http://farm8.staticflickr.com/7311/12060748443_6b1736a453_b.jpg" width="400">


TODO
======
- Check that trace drawing works, and implement exporter in OBP
- Add option setting the size of the start state node
- Add option to draw the nodes instead on edges
- Add the possibility to highlight the states in the exploration queue -- for explosion cases
	- maybe a generic approach where I get a list of nodes that I want to highlight it suffices
- Create a UI that enables seeing the contents of the states
- instead of starting with all nodes in 0,0 start them at random position 
