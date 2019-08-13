set dir "$pdb_dir"

$load_pdb_frames_into_vmd

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

