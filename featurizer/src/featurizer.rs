use std::env;
use bioshell_interactions::BackboneHBondMap;
use bioshell_io::{open_file, read_whitespace_delimited_values};
use clap::{Parser};

use bioshell_pdb::{Structure, Deposit, code_and_chain, find_cif_file_name, find_pdb_file_name, PDBError};
use bioshell_seq::chemical::{MonomerType, StandardResidueType};
use log::{debug, info, warn};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Command line application to create input data for training deep_bbq v.2 model
/// say featurizer -h to see options
struct Args {
    /// a single CIF or PDB file to process
    #[clap(short, long,  short='i')]
    input_file: Option<String>,
    /// select chain to process from the input file provided with -i option
    #[clap(short, long,  short='c')]
    select_chain: Option<String>,
    /// file with a list of PDB IDs
    #[clap(short, long,  short='l')]
    list_file: Option<String>,
    /// path to the folder with mmCIF files
    #[clap(short, long, default_value = "", short='p')]
    path: String,
}

fn find_deposit_files(list_file: &str, path: &str) -> Vec<(String, Option<String>)> {

    let reader = open_file(list_file).expect(&format!("Can't open {} file!", list_file));
    let lines: Vec<Vec<String>> = read_whitespace_delimited_values(reader).expect("Can't parse a flat text file!");
    debug!("Loading a list-file: {}", list_file);
    let mut input_files: Vec<(String, Option<String>)> = Vec::new();
    for line in lines {
        if line.len() < 1 { continue; }
        let (pdb_code, chain_id) = code_and_chain(&line[0]);
        if let Ok(cif_fname) = find_cif_file_name(&pdb_code, path) {
            input_files.push((cif_fname, chain_id));
            continue;
        }
        if let Ok(pdb_fname) = find_pdb_file_name(&pdb_code, path) {
            input_files.push((pdb_fname, chain_id));
            continue;
        }
        warn!("Can't find a PDB file for the following PDB ID: {:?}!\nSpecify folder with --path option", &pdb_code);
    }
    info!("{} input files found in {}",input_files.len(), list_file);

    return input_files;
}

fn process_deposit(fname: &str, chain: &str) -> Result<(), PDBError> {

    let deposit = Deposit::from_file(fname)?;
    let mut strctr = deposit.structure();
    strctr.remove_ligands();
    let strctr = Structure::from_iterator(&strctr.id_code, strctr.atoms().iter().filter(|a| a.chain_id == chain));
    let entity_id = &strctr.atoms()[0].entity_id;
    let entity = deposit.entity(entity_id);
    let entity_resids = entity.chain_monomers(chain)?;
    let chain_resids = strctr.residue_ids();
    let mut i_res_idx = 0;
    let hbonds = BackboneHBondMap::new(&strctr);
    for res in entity_resids {
        if res.parent_type==StandardResidueType::GAP {
            println!("{:^4} {}", '-', res);
            continue;
        }
        let i_res = &chain_resids[i_res_idx];
        if let Ok(ca) = strctr.atom(i_res, " CA ") {
            print!("{:4} {} {} : {:8.3} {:8.3} {:8.3}", i_res_idx, res, i_res, ca.pos.x, ca.pos.y, ca.pos.z);
            for (j_res_idx, j_res) in chain_resids.iter().enumerate() {
                if let Some(hb) = hbonds.get_h_bond(i_res, j_res) {
                    print!(" {:4} {:.3}", j_res_idx, hb.dssp_energy());
                }
                if let Some(hb) = hbonds.get_h_bond(j_res, i_res) {
                    print!(" {:4} {:.3}", j_res_idx, hb.dssp_energy());
                }
            }
            println!();
        } else {
            warn!("CA atom missing for residue: {}", i_res);
        }

        i_res_idx += 1;
    }
    Ok(())
}

fn main() -> Result<(), PDBError> {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    let mut input_files: Vec<(String, Option<String>)> = vec![];

    // ---------- Load a list of PDB IDs and try to locate all the files
    if let Some(fname) = args.list_file {
        input_files = find_deposit_files(&fname, &args.path);
    } else if let Some(fname) = args.input_file {
        input_files.push((fname, args.select_chain));
    } else {
        panic!("No input file provided! Use -i or -l options to specify an input file!");
    }

    for (fname, chain) in input_files {
        if let Some(chain) = chain {
            process_deposit(&fname, &chain)?;
        }
    }

    return Ok(());
}