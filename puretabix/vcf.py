from typing import Dict, Generator, Iterable, Optional, Tuple

from .fsm import (
    FSMachine,
    RegexTransition,
    SetInTransition,
    SetNotInTransition,
)


class VCFLine:
    """
    Representation of a single VCF line.

    Comment lines (starting #) will have:
      comment_raw

    Meta-information (starting ##) will have
      comment_key
      comment_value_str
      comment_value_dict

    Data lines will have:
        chrom       VCF CHROM
        pos         VCF POS
        _id         VCF IDs (if any)
        ref         VCF REF
        alt         VCF ALTs (if any)
        qual        floating point representation of VCF QUAL (if possible)
        qual_str    raw string representation of VCF QUAL
        _filter     VCF FILTERs (if any)
        info        VCF INFO as a dictionary

    Data lines may have:
        sample      Iterable of dictionaries for the samples

    For output, a VCFLine can be converted to a string.

    For convenience, several additional features exist

        is_comment  checks if either comment_raw or comment_key is defined
        get_gt()    returns the ref and/or alt alleles for each sample as a tuple of tuples of strings
    """

    comment_raw: str  # non-empty if line starts with #
    comment_key: str  # non-empty if line starts with ##
    comment_value_str: str  # non-empty if line starts with ##
    comment_value_dict: Dict[str, str]
    chrom: str
    pos: int
    _id: Tuple[str]
    ref: str
    alt: Tuple[str]
    qual: Optional[float]
    qual_str: str
    _filter: Tuple[str]
    info: Dict[str, str]
    sample: Tuple[Dict[str, str]]

    __slots__ = (
        "comment_raw",
        "comment_key",
        "comment_value_str",
        "comment_value_dict",
        "chrom",
        "pos",
        "_id",
        "ref",
        "alt",
        "qual",
        "qual_str",
        "_filter",
        "info",
        "sample",
    )

    def __init__(
        self,
        comment_raw: str,
        comment_key: str,
        comment_value_str: str,
        comment_value_dict: Dict[str, str],
        chrom: str,
        pos: int,
        _id: Iterable[str],
        ref: str,
        alt: Iterable[str],
        qual_str: str,
        _filter: Iterable[str],
        info: Dict[str, str],
        sample: Iterable[Dict[str, str]],
    ):
        # TODO validation of permitted characters
        self.comment_raw = comment_raw
        self.comment_key = comment_key
        self.comment_value_str = comment_value_str
        self.comment_value_dict = comment_value_dict
        self.chrom = chrom
        self.pos = pos
        self._id = tuple(_id)
        self.ref = ref
        self.alt = tuple(alt)
        self.qual_str = qual_str
        self._filter = tuple(_filter)
        self.info = info
        # this may be zero to many
        self.sample = tuple(sample)

        try:
            self.qual = float(qual_str)
        except ValueError:
            # if we can't do the conversion, don't worry
            # value will be None but original in qual_str
            self.qual = None

    def __str__(self) -> str:
        # CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT) (SAMPLE) (SAMPLE) ...
        if self.comment_raw:
            return "#" + self.comment_raw
        elif self.comment_key and self.comment_value_str:
            # meta information line
            # ##fileformat=VCFv4.2
            return "".join(("##", self.comment_key, "=", self.comment_value_str))
        elif self.comment_key and self.comment_value_dict:
            # meta information line
            # ##FILTER=<ID=ID,Description="description">
            return "".join(
                (
                    "##",
                    self.comment_key,
                    "=<",
                    ",".join(
                        (
                            f"{key}={value}" if value else key
                            for key, value in self.comment_value_dict.items()
                        )
                    ),
                    ">",
                )
            )
        else:
            # data line
            required = "\t".join(
                (
                    self.chrom,
                    str(self.pos),
                    ";".join(self._id),
                    self.ref,
                    ",".join(self.alt),
                    self.qual_str,  # use non-lossy stored version
                    ";".join(self._filter),
                    ";".join(
                        (
                            f"{key}={','.join(value)}" if value else key
                            for key, value in self.info.items()
                        )
                    ),
                )
            )

            if not self.sample:
                return required
            else:
                # get the superset of all keys of all samples
                keyset = set()
                keylist = []
                for sample in self.sample:
                    for key in sample.keys():
                        if key not in keyset:
                            keylist.append(key)
                            keyset.add(key)
                sample_keys = ":".join(keylist)
                line = "\t".join((required, sample_keys))

                # get the values for each key in superset or .
                for sample in self.sample:
                    values = []
                    for key in keylist:
                        values.append(sample.get(key, "."))
                    line = line + "\t" + ":".join(values)
                return line

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}({repr(self.comment_raw)},{repr(self.comment_key)},{repr(self.comment_value_str)},"
            + f"{repr(self.comment_value_dict)},{repr(self.chrom)},{repr(self.pos)},{repr(self._id)},{repr(self.ref)},{repr(self.alt)},"
            + f"{repr(self.qual_str)},{repr(self._filter)},{repr(self.info)},{repr(self.sample)})"
        )

    @property
    def is_comment(self) -> bool:
        return bool(self.comment_raw or self.comment_key)

    def get_genotype(self) -> Tuple[Tuple[str, ...], ...]:
        """
        Convenience function to get the GT of each sample, and use that to match the allels in ref & alt.

        Returns a tuple of strings.
        """
        if self.is_comment:
            raise ValueError("Cannot get_gt() of a comment")

        result = [tuple() for _ in range(len(self.sample))]
        refalt = (self.ref,) + self.alt
        for i, sample in enumerate(self.sample):
            if "GT" in sample:
                # split into separate alleles
                allele_pos = sample["GT"].replace("|", "/").split("/")
                # for each allele either keep a dot or do ref+alt lookup
                alleles = [i if i == "." else refalt[int(i)] for i in allele_pos]
                result[i] = tuple(alleles)
        return tuple(result)

    @classmethod
    def as_comment_raw(cls, comment_raw: str) -> "VCFLine":
        return cls(
            comment_raw,
            comment_key="",
            comment_value_str="",
            comment_value_dict={},
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_comment_key_string(cls, key: str, value_str: str) -> "VCFLine":
        # for example ##fileformat=VCFv4.3
        return cls(
            comment_raw="",
            comment_key=key,
            comment_value_str=value_str,
            comment_value_dict={},
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_comment_key_dict(cls, key: str, value_dict: Dict[str, str]) -> "VCFLine":
        # for example
        # ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
        return cls(
            comment_raw="",
            comment_key=key,
            comment_value_str="",
            comment_value_dict=value_dict,
            chrom="",
            pos=0,
            _id=[],
            ref="",
            alt=[],
            qual_str="",
            _filter="",
            info={},
            sample=[],
        )

    @classmethod
    def as_data(
        cls,
        chrom: str,
        pos: int,
        _id: Iterable[str],
        ref: str,
        alt: Iterable[str],
        qual_str: str,
        _filter: Iterable[str],
        info: Dict[str, str],
        sample: Iterable[Dict[str, str]],
    ):
        return cls(
            comment_raw="",
            comment_key="",
            comment_value_str="",
            comment_value_dict={},
            chrom=chrom,
            pos=pos,
            _id=_id,
            ref=ref,
            alt=alt,
            qual_str=qual_str,
            _filter=_filter,
            info=info,
            sample=sample,
        )


class VCFAccumulator:
    def __init__(self):
        self.reset()

    def reset(self):
        self.__characters = []

        self.comment_raw = ""
        self.comment_key = ""
        self.comment_value = ""
        self.comment_struct = {}
        self._comment_struct_key = ""
        self._comment_struct_value = ""
        self.chrom = ""
        self.pos = 0
        self._id = []
        self.ref = ""
        self.alt = []
        self._alt_option = ""
        self.qual = ""
        self._filter = []
        self.info = {}
        self.__infokey = ""
        self.format = []
        self.samples = []

    def to_vcfline(self):
        return VCFLine(
            self.comment_raw,
            self.comment_key,
            self.comment_value,
            self.comment_struct,
            self.chrom,
            self.pos,
            self._id,
            self.ref,
            self.alt,
            self.qual,
            self._filter,
            self.info,
            self.samples,
        )

    def append_character(self, _char):
        self.__characters.append(_char)

    def end_comment(self, _char):
        self.comment_raw = "".join(self.__characters)
        self.__characters = []

    def comment_value_to_end(self, _char):
        self.comment_value = "".join(self.__characters)
        self.__characters = []

    def comment_key_to_comment_value(self, _char):
        self.comment_key = "".join(self.__characters)
        self.__characters = []

    def comment_value_to_comment_struct_key(self, _char):
        self.comment_value = "".join(self.__characters)
        self.__characters = []

    def comment_struct_key_to_comment_struct_value(self, _char):
        comment_struct_key = "".join(self.__characters)
        assert self._comment_struct_key not in self.comment_struct, (
            self._comment_struct_key,
            self.comment_struct,
        )
        self.comment_struct[comment_struct_key] = None
        self.__characters = []

    def comment_struct_value_to_comment_struct_key(self, _char):
        self._comment_struct_value = "".join(self.__characters)
        self.__characters = []
        # this needs dict keys to be insertion ordered
        self.comment_struct[
            tuple(self.comment_struct.keys())[-1]
        ] = self._comment_struct_value

    def chrom_to_pos(self, _char):
        self.chrom = "".join(self.__characters)
        self.__characters = []

    def pos_to_id(self, _char):
        self.pos = int("".join(self.__characters))
        self.__characters = []

    def id_to_ref(self, _char):
        self._id.append("".join(self.__characters))
        self.__characters = []

    def ref_to_alt(self, _char):
        self.ref = "".join(self.__characters)
        self.__characters = []

    def alt_to_alt(self, _char):
        self.alt.append("".join(self.__characters))
        self.__characters = []

    def alt_to_qual(self, _char):
        self.alt.append("".join(self.__characters))
        self.__characters = []

    def qual_to_filter(self, _char):
        self.qual = "".join(self.__characters)
        self.__characters = []

    def filter_to_info_key(self, _char):
        self._filter.append("".join(self.__characters))
        self.__characters = []

    def info_key_to_info_key(self, _char):
        self.__infokey = "".join(self.__characters)
        self.info[self.__infokey] = None
        self.__characters = []

    def info_key_to_info_value(self, _char):
        self.__infokey = "".join(self.__characters)
        self.__characters = []

    def info_value_to_format(self, _char):
        __infovalue = "".join(self.__characters)
        if self.__infokey not in self.info:
            self.info[self.__infokey] = []
        self.info[self.__infokey].append(__infovalue)
        self.__characters = []

    def format_to_sample(self, _char):
        self.format.append("".join(self.__characters))
        self.__characters = []

    def sample_to_sample(self, _char):
        sample_str = "".join(self.__characters)
        sample_parts = sample_str.split(":")
        sample = {}
        for key, value in zip(self.format, sample_parts):
            sample[key] = value
        self.samples.append(sample)
        self.__characters = []


# states for the finite state machine to parse comments
LINE_START = "LINE_START"
COMMENT = "COMMENT"
COMMENT_KEY = "COMMENT_KEY"
COMMENT_VALUE = "COMMENT_VALUE"
COMMENT_STRUCT_KEY = "COMMENT_STRUCT_KEY"
COMMENT_STRUCT_VALUE = "COMMENT_STRUCT_VALUE"
COMMENT_STRUCT_VALUE_QUOTED = "COMMENT_STRUCT_VALUE_QUOTED"

# ##CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT) (SAMPLE)*
CHROM = "CHROM"
POS = "POS"
ID = "ID"
REF = "REF"
ALT = "ALT"
QUAL = "QUAL"
FILTER = "FILTER"
INFO_KEY = "INFO_KEY"
INFO_VALUE = "INFO_VALUE"
FORMAT = "FORMAT"
SAMPLE = "SAMPLE"


def get_vcf_fsm():

    fsm_vcf = FSMachine()

    fsm_vcf.add_transition(
        LINE_START, CHROM, SetNotInTransition, "#", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        CHROM, CHROM, SetNotInTransition, "\t", VCFAccumulator.append_character
    )
    # TODO verify CHROM in contig or assembly from comments?
    fsm_vcf.add_transition(
        CHROM, POS, SetInTransition, "\t", VCFAccumulator.chrom_to_pos
    )
    fsm_vcf.add_transition(
        POS, POS, SetInTransition, "0123456789", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(POS, ID, SetInTransition, "\t", VCFAccumulator.pos_to_id)
    fsm_vcf.add_transition(
        ID, ID, RegexTransition, r"\S", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(ID, ID, SetInTransition, ";", VCFAccumulator.id_to_ref)
    fsm_vcf.add_transition(ID, REF, SetInTransition, "\t", VCFAccumulator.id_to_ref)
    fsm_vcf.add_transition(
        REF, REF, SetInTransition, "ACGTN", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(REF, ALT, SetInTransition, "\t", VCFAccumulator.ref_to_alt)
    fsm_vcf.add_transition(
        ALT, ALT, SetNotInTransition, ",\t", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(ALT, ALT, SetInTransition, ",", VCFAccumulator.alt_to_alt)
    fsm_vcf.add_transition(ALT, QUAL, SetInTransition, "\t", VCFAccumulator.alt_to_qual)
    fsm_vcf.add_transition(
        QUAL, QUAL, SetInTransition, "0123456789.-", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        QUAL, FILTER, SetInTransition, "\t", VCFAccumulator.qual_to_filter
    )
    fsm_vcf.add_transition(
        FILTER, FILTER, SetNotInTransition, "\t;", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        FILTER, FILTER, SetInTransition, ";", VCFAccumulator.filter_to_info_key
    )
    fsm_vcf.add_transition(
        FILTER, INFO_KEY, SetInTransition, "\t", VCFAccumulator.filter_to_info_key
    )
    fsm_vcf.add_transition(
        INFO_KEY,
        INFO_KEY,
        SetNotInTransition,
        "=\t;\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        INFO_KEY,
        INFO_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.info_key_to_info_value,
    )
    fsm_vcf.add_transition(
        INFO_KEY, INFO_KEY, SetInTransition, ";", VCFAccumulator.info_key_to_info_key
    )
    fsm_vcf.add_transition(
        INFO_KEY, FORMAT, SetInTransition, "\t", VCFAccumulator.info_key_to_info_key
    )
    fsm_vcf.add_transition(
        INFO_KEY,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.info_key_to_info_key,
    )
    fsm_vcf.add_transition(
        INFO_VALUE,
        INFO_VALUE,
        SetNotInTransition,
        "\t;,\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        INFO_VALUE,
        INFO_VALUE,
        SetInTransition,
        ",",
        VCFAccumulator.info_value_to_format,
    )
    fsm_vcf.add_transition(
        INFO_VALUE, INFO_KEY, SetInTransition, ";", VCFAccumulator.info_value_to_format
    )
    fsm_vcf.add_transition(
        INFO_VALUE, FORMAT, SetInTransition, "\t", VCFAccumulator.info_value_to_format
    )
    fsm_vcf.add_transition(
        INFO_VALUE,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.info_value_to_format,
    )
    fsm_vcf.add_transition(
        FORMAT, FORMAT, SetNotInTransition, "\t:", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        FORMAT, FORMAT, SetInTransition, ":", VCFAccumulator.format_to_sample
    )
    fsm_vcf.add_transition(
        FORMAT, SAMPLE, SetInTransition, "\t", VCFAccumulator.format_to_sample
    )
    fsm_vcf.add_transition(
        SAMPLE, SAMPLE, SetNotInTransition, "\t\n", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        SAMPLE, SAMPLE, SetInTransition, "\t", VCFAccumulator.sample_to_sample
    )
    fsm_vcf.add_transition(
        SAMPLE, None, SetInTransition, ("\n", None), VCFAccumulator.sample_to_sample
    )
    fsm_vcf.add_transition(LINE_START, COMMENT, SetInTransition, "#", None)
    fsm_vcf.add_transition(
        COMMENT, COMMENT, SetNotInTransition, "#\n", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        COMMENT, None, SetInTransition, ("\n", None), VCFAccumulator.end_comment
    )
    fsm_vcf.add_transition(COMMENT, COMMENT_KEY, SetInTransition, "#", None)
    fsm_vcf.add_transition(
        COMMENT_KEY,
        COMMENT_KEY,
        SetNotInTransition,
        "=",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_KEY,
        COMMENT_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.comment_key_to_comment_value,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        COMMENT_VALUE,
        SetNotInTransition,
        "<\n",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        COMMENT_STRUCT_KEY,
        SetInTransition,
        "<",
        VCFAccumulator.comment_value_to_comment_struct_key,
    )
    fsm_vcf.add_transition(
        COMMENT_VALUE,
        None,
        SetInTransition,
        ("\n", None),
        VCFAccumulator.comment_value_to_end,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_KEY,
        COMMENT_STRUCT_KEY,
        SetNotInTransition,
        "=",
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_KEY,
        COMMENT_STRUCT_VALUE,
        SetInTransition,
        "=",
        VCFAccumulator.comment_struct_key_to_comment_struct_value,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_VALUE,
        SetNotInTransition,
        '",>',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_VALUE_QUOTED,
        SetInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE_QUOTED,
        COMMENT_STRUCT_VALUE_QUOTED,
        SetNotInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE_QUOTED,
        COMMENT_STRUCT_VALUE,
        SetInTransition,
        '"',
        VCFAccumulator.append_character,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT_STRUCT_KEY,
        SetInTransition,
        ",",
        VCFAccumulator.comment_struct_value_to_comment_struct_key,
    )
    fsm_vcf.add_transition(
        COMMENT_STRUCT_VALUE,
        COMMENT,
        SetInTransition,
        ">",
        VCFAccumulator.comment_struct_value_to_comment_struct_key,
    )

    return fsm_vcf


def read_vcf_lines(input: Iterable[str]) -> Generator[VCFLine, None, None]:
    """
    Convenience function for parsing a source of VCF lines
    """
    vcf_fsm = get_vcf_fsm()
    accumulator = VCFAccumulator()
    for line in input:
        vcf_fsm.run(line, LINE_START, accumulator)
        vcfline = accumulator.to_vcfline()
        accumulator.reset()
        yield vcfline
