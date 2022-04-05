use dot_writer::{Attributes, Color, DotWriter, Shape, Style};
use regex::RegexSetBuilder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Hash, Eq, PartialEq, Debug, Clone)]
struct Node {
    name: String,
    cons: Vec<String>,
}

fn main() {
    let log_path = "/home/brad/gigabit_code/log/output.log";
    let buffered = BufReader::new(File::open(log_path).unwrap());

    let set = RegexSetBuilder::new(&[r#"TRACE"#])
        .case_insensitive(true)
        .build()
        .unwrap();
    let mut trace_lines = Vec::new();

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
        for ind in 2..words.len() {
            list.push(words[ind].to_string());
        }
    }
    // println!("{:?}", list);

    let mut funcs = Vec::<String>::new();
    let mut verbs = Vec::<String>::new();
    let mut in_loop = false;

    for ent in list.iter() {
        if ent == "END" {
            in_loop = false;
        }
        if in_loop {
            if ((ent == "start") ^ (ent == "end")) {
                verbs.push(ent.to_string());
            } else {
                funcs.push(ent.to_string());
            }
        }
        if ent == "START" {
            in_loop = true;
        }
    }
    for ind in 0..funcs.len() {
        println!("{} {}", funcs[ind], verbs[ind]);
    }

    // nodes.insert("main", Vec::<String>::new());
    // println!("{:?}", nodes);

    let mut last_func: &String = &"".to_string();
    let mut last_verb: &String = &"".to_string();
    let mut verb = &"".to_string();
    let mut connections = Vec::new();
    let mut cluster_name = &"".to_string();
    let mut cluster_nodes: Vec<String> = Vec::new();
    let mut cluster_name2 = &"".to_string();
    let mut cluster_nodes2: Vec<String> = Vec::new();
    let mut cluster_name3 = &"".to_string();
    let mut cluster_nodes3: Vec<String> = Vec::new();
    let mut clusters_parent = HashMap::<String, Vec<String>>::new();
    let mut clusters_child = HashMap::<String, Vec<String>>::new();
    let mut clusters_childschild = HashMap::<String, Vec<String>>::new();
    let mut cluster_flag1 = false;
    let mut cluster_flag2 = false;
    let mut cluster_flag3 = false;
    let mut cluster_entry_func = &"".to_string();
    let mut child_parent = Vec::new();
    let mut childs_child = Vec::new();
    let mut new_clust_name = String::new();
    // let mut node_renames = Vec::new();

    // connections.push(["".to_string(), "".to_string()]);
    // println!("len {}, inside len {}", connections.len(), connections[0].len());

    for (ind, fcn) in funcs.iter().enumerate() {
        verb = &verbs[ind];
        // println!("{} {} {} {}, ind {}, 1 {}, 2 {}, 3 {} \n", fcn, verb,last_func, last_verb, ind, verb == &"start".to_string() &&  last_verb == &"end".to_string(), verb == &"start".to_string() &&  last_verb == &"start".to_string(), verb == &"end".to_string() &&  last_verb == &"end".to_string());

        if verb == &"start".to_string() && last_verb == &"end".to_string() {
            if cluster_flag1 == false && cluster_flag2 == false {
                // println!("connection {} {} ", last_func, fcn);
                connections.push([last_func.to_string(), fcn.to_string()]);
                last_func = fcn;
                last_verb = verb;
            } else if cluster_flag1 == true {
                if cluster_flag2 == false {
                    // connections.push([last_func.to_string(), fcn.to_string()]);
                    // println!("connection {} {} ", last_func, fcn);
                    cluster_nodes.push(fcn.to_string());
                    last_func = fcn;
                    last_verb = verb;
                } else if cluster_flag2 == true {
                    if cluster_flag3 == false {
                        // connections.push([last_func.to_string(), fcn.to_string()]);
                        // println!("connection {} {} ", last_func, fcn);
                        cluster_nodes2.push(fcn.to_string());
                        last_func = fcn;
                        last_verb = verb;
                    } else if cluster_flag3 == true {
                        // connections.push([last_func.to_string(), fcn.to_string()]);
                        last_func = fcn;
                        last_verb = verb;
                        // cluster_nodes.push(fcn.to_string());
                        cluster_nodes3.push(fcn.to_string());
                    }
                    // // connections.push([last_func.to_string(), fcn.to_string()]);
                    // last_func = fcn;
                    // last_verb = verb;
                    // // cluster_nodes.push(fcn.to_string());
                    // cluster_nodes2.push(fcn.to_string());
                }
            }
        } else if verb == &"start".to_string() && last_verb == &"start".to_string() {
            if cluster_flag1 == false {
                cluster_flag1 = true;
                cluster_name = last_func;
                cluster_nodes.push(fcn.to_string());
                cluster_entry_func = last_func;
            } else if cluster_flag1 == true && cluster_flag2 == false {
                cluster_flag2 = true;
                cluster_name2 = last_func;
                // cluster_nodes.push(fcn.to_string());
                cluster_nodes2.push(fcn.to_string());
                cluster_entry_func = last_func;
                child_parent.push((cluster_name, cluster_name2));
            } else {
                cluster_flag3 = true;
                cluster_name3 = last_func;
                // cluster_nodes.push(fcn.to_string());
                cluster_nodes3.push(fcn.to_string());
                cluster_entry_func = last_func;
                childs_child.push((cluster_name2, cluster_name3));
            }

            // connections.push([last_func.to_string(), fcn.to_string()]);

            last_func = fcn;
            last_verb = verb;
        } else if verb == &"end".to_string() && last_verb == &"end".to_string() {
            if cluster_flag3 == false {
                if cluster_flag2 == false {
                    cluster_flag1 = false;
                    clusters_parent.insert(cluster_name.to_string(), cluster_nodes.clone());
                    cluster_nodes.clear();
                } else {
                    cluster_flag2 = false;
                    clusters_child.insert(cluster_name2.to_string(), cluster_nodes2.clone());
                    cluster_nodes2.clear();
                }
            } else {
                cluster_flag3 = false;
                clusters_childschild.insert(cluster_name3.to_string(), cluster_nodes3.clone());
                cluster_nodes3.clear();
            }

            // connections.push([last_func.to_string(), fcn.to_string()]);
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
    // println!("\n\n\n parent clusters: {:?}", clusters_parent);
    // println!("\n\n\n child clusters: {:?}", clusters_child);
    // }
    println!("\nConnection \n");
    for ent in connections.iter() {
        println!("{} -> {}", ent[0], ent[1]);
    }

    println!("\nParent \n");
    for (ent, val) in clusters_parent.iter() {
        println!("{}\n", ent);
        for item in val.iter() {
            println!("\t {}", item);
        }
    }

    println!("\nChildren \n");
    for (ent, val) in clusters_child.iter() {
        println!("{}\n", ent);
        for item in val.iter() {
            println!("\t {}", item);
        }
    }

    let mut f = File::create("example3.dot").unwrap();
    let mut clusters_parent_id = HashMap::<String, Vec<dot_writer::NodeId>>::new();
    let mut clusters_child_id = HashMap::<String, Vec<dot_writer::NodeId>>::new();
    let mut clusters_childschild_id = HashMap::<String, Vec<dot_writer::NodeId>>::new();

    let mut output_bytes = Vec::new();
    {
        let mut writer = DotWriter::from(&mut output_bytes);
        writer.set_pretty_print(true);
        let mut digraph = writer.digraph();
        // digraph.set_rank_direction(dot_writer::RankDirection::LeftRight);
        {
            for (clust, nodes) in &clusters_parent {
                let mut clust_write = digraph.cluster();
                clust_write.set_label(clust);
                clust_write.set_rank_direction(dot_writer::RankDirection::LeftRight);
                let mut node_ids = Vec::new();
                for item in nodes.clone().iter() {
                    let mut nod = clust_write.node_auto();
                    nod.set_label(item);
                    node_ids.push(nod.id());
                }
                clusters_parent_id.insert(clust.to_string(), node_ids.clone());

                let length = nodes.len();
                // clust_write.edge(&clust, &node_ids[0]);

                if nodes.len() > 1 {
                    for ind in 1..nodes.len() {
                        clust_write.edge(&node_ids[ind - 1], &node_ids[ind]);
                    }
                }
                for pair in &child_parent {
                    if pair.0 == clust {
                        for (sub_clust, sub_nodes) in &clusters_child {
                            if sub_clust == pair.1 {
                                let sub_clust_name = pair.1;
                                let mut sub_clust = clust_write.cluster();
                                sub_clust.set_label(sub_clust_name);
                                let mut node_ids2 = Vec::new();
                                for item in sub_nodes.clone().iter() {
                                    let mut nod2 = sub_clust.node_auto();
                                    nod2.set_label(item);
                                    node_ids2.push(nod2.id());
                                }
                                clusters_child_id
                                    .insert(sub_clust_name.to_string(), node_ids2.clone());
                                let length2 = node_ids2.len();
                                // sub_clust.edge(&sub_clust_name, &node_ids2[0]);
                                for ind2 in 1..length2 {
                                    sub_clust.edge(&node_ids2[ind2 - 1], &node_ids2[ind2]);
                                }
                                for pair2 in &childs_child {
                                    if pair2.0 == sub_clust_name {
                                        for (subsub_clust, subsub_nodes) in &clusters_childschild {
                                            if subsub_clust == pair2.1 {
                                                let subsub_clust_name = pair2.1;
                                                let mut subsub_clust = sub_clust.cluster();
                                                subsub_clust.set_label(subsub_clust_name);
                                                let mut node_ids3 = Vec::new();
                                                for item2 in subsub_nodes.clone().iter() {
                                                    let mut nod3 = subsub_clust.node_auto();
                                                    nod3.set_label(item2);
                                                    node_ids3.push(nod3.id());
                                                }
                                                clusters_childschild_id.insert(
                                                    subsub_clust_name.to_string(),
                                                    node_ids3.clone(),
                                                );

                                                let length3 = node_ids3.len();
                                                // sub_clust.edge(&sub_clust_name, &node_ids2[0]);
                                                for ind3 in 1..length3 {
                                                    subsub_clust.edge(
                                                        &node_ids3[ind3 - 1],
                                                        &node_ids3[ind3],
                                                    );
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            for item in connections.iter() {
                digraph.edge(item[0].to_string(), item[1].to_string());
            }
        }
        // digraph goes out of scope here to write closing bracket
        // then writer goes out of scope here to free up output_bytes for reading
    }

    println!(
        "echo \"{}\" | dot > ex.dot",
        String::from_utf8(output_bytes).unwrap()
    );
}
