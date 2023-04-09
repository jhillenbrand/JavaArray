package net.sytes.botg.datatypes;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

@Retention(RetentionPolicy.RUNTIME) // retention policy must be set to runtime, in order for parser to find it while configuration
@Target(ElementType.FIELD)
public @interface ConfigOption {
	String description() default "";
}

