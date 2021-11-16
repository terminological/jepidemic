package org.github.terminological.jepidemic.growth;

import java.util.HashMap;
import java.util.Optional;

public class OptionalMap<K, V> extends HashMap<K, V> {

	public Optional<V> optional(K key) {return Optional.ofNullable(get(key));}
	
}
