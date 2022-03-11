use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::RegexSetBuilder;
use std::collections::HashMap;
use dot_writer::{Color, DotWriter, Attributes, Shape, Style};


#[derive(Hash, Eq, PartialEq, Debug, Clone)]
struct Node{
    name: String,
    cons: Vec<String>,
}

fn main() {
    let log_path = "/home/brad/gigabit_code/log/output.log";
    let buffered = BufReader::new(File::open(log_path).unwrap());

    let set = RegexSetBuilder::new(&[
        r#"TRACE"#,
    ]).case_insensitive(true)
        .build().unwrap();
    let mut trace_lines= Vec::new();

    buffered
        .lines()
        .filter_map(|line| line.ok())
        .filter(|line| set.is_match(line.as_str()))
        // .for_each(|x| println!("{}", x));
        .for_each(|x| trace_lines.push(x));


    // Splits all trace entries from log into a vector of strings containing the information after TRACE -, does so one word at a time 
    let mut list: Vec<String> = Vec::new();

    for x in trace_lines.iter() {
        // let childvec = Vec::new();
        let words: Vec<&str> = x.trim().split(" ").collect();
        for ind in 2 .. words.len(){
            list.push(words[ind].to_string());
        }  
    }
    println!(" list {:?}", list);

    let mut funcs = Vec::<String>::new();
    let mut verbs = Vec::<String>::new();
    let mut in_loop = false;

    for ent in list.iter(){

        if ent == "END"{
            in_loop = false;
        }
        if in_loop{
            if ((ent == "start") ^ (ent == "end")){
                verbs.push(ent.to_string());
            } else{
                funcs.push(ent.to_string());
            } 
        } 
        if ent == "START"{
            in_loop = true;
        }           
    }
    for ind in 0 .. funcs.len(){
        println!("{} {}", funcs[ind], verbs[ind]);
    }
    

    // nodes.insert("main", Vec::<String>::new());
    // println!("{:?}", nodes);
    let mut nodes: Vec<Vec<&String>> = Vec::new();
    let mut last_func: &String = &"".to_string();
    let mut last_verb: &String = &"".to_string();
    let mut verb = &"".to_string();
    let mut connections = Vec::new();
    let mut cluster_name = &"".to_string();
    let mut cluster_nodes: Vec<String> = Vec::new();
    let mut cluster_name2 = &"".to_string();
    let mut cluster_nodes2: Vec<String> = Vec::new();
    let mut clusters_parent = HashMap::<String, Vec<String>>::new();
    let mut clusters_child = HashMap::<String, Vec<String>>::new();
    let mut cluster_flag1 = false;
    let mut cluster_flag2 = false;
    let mut cluster_entry_func = &"".to_string();
    let mut child_parent = Vec::new();
    let mut node_id = String::new();


    // connections.push(["".to_string(), "".to_string()]);
    // println!("len {}, inside len {}", connections.len(), connections[0].len());
        
    for (ind,fcn) in funcs.iter().enumerate(){
        verb = &verbs[ind];

        // println!("{} {} {} {}, ind {}, 1 {}, 2 {}, 3 {} \n", fcn, verb,last_func, last_verb, ind, verb == &"start".to_string() &&  last_verb == &"end".to_string(), verb == &"start".to_string() &&  last_verb == &"start".to_string(), verb == &"end".to_string() &&  last_verb == &"end".to_string());    

        if verb == &"start".to_string() &&  last_verb == &"end".to_string(){
            // println!("connection {} {} ", last_func, fcn);

            connections.push([last_func.to_string(), fcn.to_string()]);
            if cluster_flag2 {
                cluster_nodes2.push(node_id);
            } 
            if cluster_flag1 {
                cluster_nodes.push(node_id);
            }
            last_func = fcn;
            last_verb = verb;

        } else if verb == &"start".to_string() &&  last_verb == &"start".to_string(){
            if cluster_flag1 == false{
                cluster_flag1 = true;
                cluster_name = last_func;
                cluster_nodes.push(node_id);
                cluster_entry_func = last_func;
            } else {
                cluster_flag2 = true;
                cluster_name2 = last_func;
                cluster_nodes.push(node_id);
                cluster_nodes2.push(node_id);
                cluster_entry_func = last_func;
                child_parent.push((cluster_name, cluster_name2));
            }
            
            // connections.push([last_func.to_string(), fcn.to_string()]);
            
            last_func = fcn;
            last_verb = verb;
        } else if verb == &"end".to_string() &&  last_verb == &"end".to_string(){
            if cluster_flag2 == false{
                cluster_flag1 = false;
                clusters_parent.insert(cluster_name.to_string(), cluster_nodes.clone());
                cluster_nodes.clear();
            } else {
                cluster_flag2 = false;
                clusters_child.insert(cluster_name2.to_string(), cluster_nodes2.clone());
                cluster_nodes2.clear();
            }
            connections.push([last_func.to_string(), node_id]);
            last_func = fcn;
            last_verb = verb;
            
        } else {
            last_func = fcn;
            last_verb = verb;
        }

        //set last func/verb
        

    }   

    // for item in connections.clone().iter(){
    //     println!("\n {} -> {}", item[0], item[1]);
    // }

    // sort and remove duplicates
    let mut functions = funcs.clone();
    // println!("funcs: {:?}", functions);
    functions.sort();
    functions.dedup();
    // println!("funcs: {:?}", functions);
    // {
    // println!("connections: {:?}", connections);
    connections.sort();
    connections.dedup();
    // println!("connections: {:?}", connections);
    // }

    let mut f = File::create("example3.dot").unwrap();

    let mut output_bytes = Vec::new();
{
    let mut writer = DotWriter::from(&mut output_bytes);
    writer.set_pretty_print(true);
    let mut digraph = writer.digraph();
    {
        for (clust, nodes) in &clusters_parent{
            let mut clust_write = digraph.cluster();
            clust_write.set_label(clust);
            let length = nodes.len();
            clust_write.edge(&clust, &nodes[length-1]);
            if nodes.len() > 1{
                for ind in 1 .. nodes.len(){
                    clust_write.edge(&nodes[ind-1], &nodes[ind]);
                }
            }
            for pair in &child_parent{
                if pair.0 == clust{
                    for (sub_clust, sub_nodes) in &clusters_child{
                        if sub_clust == pair.1{
                            let sub_clust_name = pair.1;
                            let mut sub_clust = clust_write.subgraph();
                            let length2 = sub_nodes.len();
                            sub_clust.edge(&sub_clust_name, &sub_nodes[0]);
                            for ind2 in 1 .. length2 {
                                sub_clust.edge(&sub_nodes[ind2-1], &sub_nodes[ind2]);
                            }
                        }
                    }
                    
                }
            }
        }

        for item in connections.iter(){
            digraph.edge(item[0].to_string(), item[1].to_string());
        }

    }
    // digraph goes out of scope here to write closing bracket
    // then writer goes out of scope here to free up output_bytes for reading
}

println!("echo \"{}\" | dot > ex.dot", String::from_utf8(output_bytes).unwrap());

let mut output_bytes2 = Vec::new();
{
let mut writer2 = DotWriter::from(&mut output_bytes2);
    writer2.set_pretty_print(true);
    let mut digraph = writer2.digraph();
    digraph.node_named("a");
    {
        let mut clust5 = digraph.cluster();
        clust5.set_label("B");
        clust5.set_rank_direction(dot_writer::RankDirection::TopBottom);
        clust5.node_named("b").set_rank(dot_writer::Rank::Min);
        clust5.node_named("c").set_rank(dot_writer::Rank::Max);
        clust5.node_named("d").set_rank(dot_writer::Rank::Max);
        clust5.edge("c", "d");
        {
            let mut clust7 = clust5.cluster();
            clust7.cluster();
            clust7.set_label("E");
            clust7.node_named("e");
            clust7.node_named("f");
            clust7.edge("e", "f"); 
        }
        clust5.node_named("i");
        {
            let mut clust7 = clust5.cluster();
            clust7.cluster();
            clust7.set_label("E");
            clust7.node_named("j");
            clust7.node_named("k");
            clust7.edge("f", "j"); 
        }
    }
digraph.edge("a","b");
digraph.edge("b","g");
digraph.edge("b", "e");
}
    println!("echo \"{}\" | dot > ex.dot", String::from_utf8(output_bytes2).unwrap());
}
    
impl Node{
    fn new(ent: &str) -> Node{
        let node = Node{
            name: ent.to_string(),
            cons: Vec::<String>::new(),
        };
        node
    }
}