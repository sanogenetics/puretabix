from .fsm import (
    FSMachine,
    RegexTransition,
    SetInTransition,
    SetNotInTransition,
)


class VCFLine:

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
        comment_value_dict,
        chrom: str,
        pos: int,
        _id,
        ref: str,
        alt,
        qual_str: str,
        _filter,
        info,
        sample,
    ):
        self.comment_raw = comment_raw
        self.comment_key = comment_key
        self.comment_value_str = comment_value_str
        self.comment_value_dict = comment_value_dict
        self.chrom = chrom
        self.pos = pos
        self._id = _id
        self.ref = ref
        self.alt = alt
        self.qual_str = qual_str
        self._filter = _filter
        self.info = info
        # this may be zero to many
        self.sample = sample

        try:
            self.qual = float(qual_str)
        except ValueError:
            # if we can't do the conversion, don't worry
            # value will be None but original in qual_str
            self.qual = None

    def __str__(self):
        # CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT) (SAMPLE) (SAMPLE) ...
        if self.comment_raw:
            return "#" + self.comment_raw
        elif self.comment_key and self.comment_value_str:
            return "".join(("##", self.comment_key, "=", self.comment_value_str))
        elif self.comment_key and self.comment_value_dict:
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
                sample_keys = ":".join(self.sample[0].keys())
                line = "\t".join((required, sample_keys))

                for sample in self.sample:
                    sample_values = ":".join(sample.values())
                    line = line + "\t" + sample_values
                return line

    @classmethod
    def as_comment_raw(cls, comment_raw):
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
    def as_comment_key_string(cls, key, value_str):
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
    def as_comment_key_dict(cls, key, value_dict):
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
        ID, ID, RegexTransition, r"[a-zA-Z0-9.]", VCFAccumulator.append_character
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
        INFO_KEY, None, SetInTransition, "\n", VCFAccumulator.info_key_to_info_key
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
        INFO_VALUE, None, SetInTransition, "\n", VCFAccumulator.info_value_to_format
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
        SAMPLE, None, SetInTransition, "\n", VCFAccumulator.sample_to_sample
    )
    fsm_vcf.add_transition(LINE_START, COMMENT, SetInTransition, "#", None)
    fsm_vcf.add_transition(
        COMMENT, COMMENT, SetNotInTransition, "#\n", VCFAccumulator.append_character
    )
    fsm_vcf.add_transition(
        COMMENT, None, SetInTransition, "\n", VCFAccumulator.end_comment
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
        COMMENT_VALUE, None, RegexTransition, "\n", VCFAccumulator.comment_value_to_end
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
