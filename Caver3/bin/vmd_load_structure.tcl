#set dir "$pdb_dir"

mol load pdb ../data/$pdb_representant

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

