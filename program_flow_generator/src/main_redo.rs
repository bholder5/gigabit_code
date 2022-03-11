use std::fs::File;
use std::io::{BufReader, BufRead};
use regex::RegexSetBuilder;
use std::collections::HashMap;
use dot_writer::{Color, DotWriter, Attributes, Shape, Style};

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
    // println!(" list {:?}", list);

    let mut funcs = Vec::<String>::new();
    let mut verbs = Vec::<String>::new();
    let mut in_loop = false;

    for ent in list.iter(){
        if ent == "END"{
            in_loop = false;
        }
        if in_loop{
            if (ent == "start") ^ (ent == "end"){
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

    let mut output_bytes = Vec::new();
    {
        let mut writer = DotWriter::from(&mut output_bytes);
        writer.set_pretty_print(true);
        let mut digraph = writer.digraph();

        let mut last_fnc: &String = &"".to_string();
        let mut last_verb: &String = &"".to_string();
        let mut verb = &"".to_string();
        let mut fnc = &"".to_string();

        for ind in 1 .. (&funcs.len()).clone(){
            last_fnc = &funcs[ind-1];
            last_verb = &verbs[ind-1];
            fnc = &funcs[ind];
            verb = &verbs[ind];

            
        }
    }      
}   