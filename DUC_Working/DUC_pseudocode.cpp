void perform_DUC(vector<cell* > initial_domain, int discretization) {

	set<cell* > current_domain; // contains all the nodes currently in the domain
	queue<cell* > next_generation_divide; // D_k in the multilevel sampling algorithm
	queue<cell* > current_flagged_divide; // all the nodes that were flagged by divide
	auto itr = initial_domain.begin();
	auto en = initial_domain.end();
	while (itr != en) {
		current_domain.insert(*itr);
		current_flagged_divide.push(*itr);
		itr++;
	}

	queue< pair<cell*,cell* > > current_flagged_united; // all the nodes that were flagged by united
	int cutoff = 0;
	bool done_divide = false, done_unite = false, done = false;
	while (!current_flagged_divide.empty()) {
		pop node from current_flagged_divide;
		remove node from current_domain;
		insert children of node into current_domain; (from the tree)
		insert children of node into next_generation_divide;
	}
	for (int i = 0; i < discretization && !done; i++) {
		// perform one iteration of the multilevel sampling,flag some nodes by populating current_flagged_divide
		if (!done_divide)
			flag_divide(next_generation_divide, current_flagged_divide, cutoff, done_divide); 
		// perform the actual division
		while (!current_flagged_divide.empty()) {
			pop node from current_flagged_divide;
			remove node from current_domain;
			insert children of node into current_domain; (from the tree)
			insert children of node into next_generation_divide;
		}
		if (!done_unite)
			flag_unite(current_domain, current_flagged_unite); 	// perform the actual unite
		while (!current_flagged_unite.empty()) {
			pair<cell*, cell* > cells_to_unite = pop current_flagged_unite;
			remove both the nodes from current_domain;
			splice these nodes from their parents;
			cell* new_parent = find_LCA(cells_to_unite);
			cell* new_node = fuse(cells_to_unite);
			delete cells_to_unite.first;
			delete cells_to_unite.second;
			new_node->parent = new_parent;

		}
		done = done_unite && done_divide;	
	}
}
