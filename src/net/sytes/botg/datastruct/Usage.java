package net.sytes.botg.datastruct;

import static java.lang.annotation.RetentionPolicy.RUNTIME;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.Target;

/**
 * Annotation that can be used to describe the purpose of a class and explain its usage
 * by providing a description text
 * @author hillenbrand
 *
 */
@Retention(RUNTIME)
@Target(ElementType.TYPE)
public @interface Usage {
	String description() default "";
}
