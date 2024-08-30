//Credit: https://www.algotree.org/algorithms/backtracking/partition_n_elements_into_k_nonempty_subsets/
#include <stdlib.h>
#include<iostream>
#include<vector>
#include <OpenMesh/Apps/Assembly/subset_division.hh>

using namespace std;



// Function for partitioning n elements into k non-empty subsets.
void Generate_Subsets (int element, int n, int k, int non_empty_subsets, 
                       vector<vector<int>>& subsets, vector<vector<vector<int> > >& all_subsets, int &subset_count) {
    
    if (element > n) { // Max available number in a set is n.
        if (non_empty_subsets == k) { // We have generted k non-empty subsets.
            //cout << "---------------" << endl;
            subset_count++;
            all_subsets.push_back(subsets);
            /*cout << "Subset : " << subset_count << endl;
            for (auto& set : subsets) {
                for(auto& num : set) {
                    cout << num << " ";
                } cout << endl;
            }*/
        }
        return;
    }

    for (int i=0; i<subsets.size(); i++) {
        if (subsets[i].size() > 0) {
            subsets[i].push_back(element);
            Generate_Subsets(element+1, n, k, non_empty_subsets, subsets, all_subsets, subset_count);
            subsets[i].pop_back();
        } else { // empty subset
            subsets[i].push_back(element);
            Generate_Subsets(element+1, n, k, non_empty_subsets+1, subsets, all_subsets, subset_count);
            subsets[i].pop_back();
            break;
        }
    }
}

//generate all possible ways to distribute 0,1,2,...,n-1 into k baskets so that no basket is left empty.
vector<vector<vector<int> > > get_partitions(int n, int k){
    vector<vector<int>> subsets(k);
    vector<vector<vector<int> > > all_subsets; 
    int subset_count = 0;

    Generate_Subsets(0, n-1, k, 0, subsets, all_subsets, subset_count); // Start by pushing the 1'st element

    return all_subsets;
}


/*
int main() {
    int n = 5; // Given N elements. ( 1..n )
    int k = 3; // To be partitioned into k subsets.

    vector<vector<vector<int> > > all_subsets = get_partitions(n, k);

    cout<<"Number of ways: "<<all_subsets.size()<<std::endl;

    for (auto subset: all_subsets){
        std::cout<<"-------------"<<std::endl;
        for (auto group: subset){
            for (auto element: group){
                std::cout<<" "<<element;
            }
            std::cout<<std::endl;
        }
        
    }

    k = 2;
    subset_count = 0;
    vector<vector<int>> subsets_1(k);   
    cout << "===================================" << endl;
    cout << "N : " << n << " K: " << k << endl;
    Generate_Subsets(1, n, k, 0, subsets_1); // Start by pushing the 1'st element
   
    return 0;
}
*/