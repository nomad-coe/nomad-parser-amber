package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object AmberParserTests extends Specification {
  "AmberParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(AmberParser, "parsers/amber/test/examples/mdout", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(AmberParser, "parsers/amber/test/examples/mdout", "json") must_== ParseResult.ParseSuccess
    }
  }
}
