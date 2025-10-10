if not IsBound( BOTT_SO ) then
    BindGlobal( "BOTT_SO" , "/home/rlutowsk/spinc-bott-fast-populate-orbit/gap.so" );
fi;

if not IsBound( BottInit ) then
    LoadDynamicModule( BOTT_SO );
fi;

BottCalcInfo := NewInfoClass("BottCalcInfo");
BottCalcInfoPrefix := "";

# Read( "/home/rlutowsk/spinc-bott-fast-populate-orbit/code.g" );
# Read( "/home/rlutowsk/spinc-bott-fast-populate-orbit/redis.g" );
Read( "/home/rlutowsk/spinc-bott-fast-populate-orbit/avl.g" );
Read( "/home/rlutowsk/spinc-bott-fast-populate-orbit/dict.g" );
